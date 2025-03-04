#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.genome = null
params.fasta = getGenomeAttribute('fasta')
params.gtf = getGenomeAttribute('gtf')
params.genome_name = null
params.genome_dir = null
params.chemistry = 'v3'
params.samp_list = null
params.concatenate = false
params.fastq = null
params.trim = false
params.skip_fastqc = false

log.info """\
    P A R S E  N E X T F L O W - I N T E G R A T E D  P I P E L I N E
    ( P A R S E N I P)
    ==================================================================
    Output folder                   : ${params.outdir}
    FASTA file(s)                   : ${params.fasta}
    GTF file(s)                     : ${params.gtf}
    Genome Name(s)                  : ${params.genome_name}
    Genome                          : ${params.genome}
    Pre-built genome index folder   : ${params.genome_dir}
    FASTQ files folder              : ${params.fastq}
    Chemistry version               : ${params.chemistry}
    Sample list file                : ${params.samp_list}
    Skip FastQC                     : ${params.skip_fastqc}
    Hard trim FASTQ files           : ${params.trim}
    Concatenate FASTQ files         : ${params.concatenate}
    Author                          : ${workflow.manifest.author}
    Homepage                        : ${workflow.manifest.homePage}
    Pipeline Version                : ${workflow.manifest.version}
    """
    .stripIndent()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHECK_SETUP_BUILD_INDEX } from './modules/local/python3_check_setup/check_setup_build_index/main.nf'
include { CHECK_SETUP_PREBUILT_INDEX } from './modules/local/python3_check_setup/check_setup_prebuilt_index/main.nf'
include { CHECK_SETUP_DUMMY_RUN } from './modules/local/python3_check_setup/check_setup_dummy_run/main.nf'

include { SPLITPIPE_INDEX } from './modules/local/splitpipe/index/main.nf'
include { SPLITPIPE_MAP } from './modules/local/splitpipe/map/main.nf'
include { CONCATENATE_FILES } from './modules/local/cat/contatenate_files/main.nf'
include { SPLITPIPE_COMBINE } from './modules/local/splitpipe/combine/main.nf'
include { TRIMGALORE_TRIM } from './modules/local/trim_galore/trim_galore_hardtrim/main.nf'
include { FASTQC } from './modules/local/fastqc/fastqc/main.nf'
include { SPLITPIPE_QC} from './modules/local/splitpipe_qc/splitpipe_qc/main.nf'
include { SCANPY_SINGLECELLQC } from './modules/local/scanpy/singlecellqc/main.nf'
include { SCANPY_SINGLECELLQC as SCANPY_SINGLECELLQC2} from './modules/local/scanpy/singlecellqc/main.nf'
include { SEURAT_BUILDSEURATOBJECT} from './modules/local/seurat/build_object/main.nf'
include { SEURAT_BUILDSEURATOBJECT as SEURAT_BUILDSEURATOBJECT2} from './modules/local/seurat/build_object/main.nf'
include { MULTIQC } from './modules/local/multiqc/multiqc/main.nf'


// Variables set in functions called below
build_index = false
perform_mapping = false


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    checkParameters()


    // TODO
    // Check input file exit and samplesheet is okay before staring the pipeline
    if(perform_mapping) {   // Only check if performing the mapping
        if(build_index) {

            fasta_ch = Channel.fromPath(params.fasta.tokenize(',')).collect()
            gtf_ch = Channel.fromPath(params.gtf.tokenize(',')).collect()

            //CHECK_SETUP_BUILD_INDEX(params.fastq, file(params.fastq), params.samp_list, file(params.samp_list), params.chemistry, params.fasta, fasta_ch, params.gtf, gtf_ch, params.genome_name)
            CHECK_SETUP_BUILD_INDEX(file(params.fastq), file(params.samp_list), params.chemistry, fasta_ch, gtf_ch, params.genome_name)
            
            check_setup_ch = CHECK_SETUP_BUILD_INDEX.out
        } else {
            CHECK_SETUP_PREBUILT_INDEX(file(params.fastq), file(params.samp_list), params.chemistry, file(params.genome_dir))
            check_setup_ch = CHECK_SETUP_PREBUILT_INDEX.out
        }
    } else {    //  Get channel output so later processes can proceed (don't check if only building a genome index)
            CHECK_SETUP_DUMMY_RUN()
            check_setup_ch = CHECK_SETUP_DUMMY_RUN.out
    }


  


    // Create split-pipe index file
    if(build_index) {
        println("Building genome index")
        
        //fasta_ch = Channel.fromPath(params.fasta.tokenize(',')).collect()
        //gtf_ch = Channel.fromPath(params.gtf.tokenize(',')).collect()
        def genome_names = params.genome_name.replace(",", " ")    // Convert comma-separated to space-separated, for split-pipe

        SPLITPIPE_INDEX(fasta_ch, gtf_ch, genome_names, check_setup_ch)
    }


     def fastq_files = [params.fastq + '/*r_{1,2}.fq.gz', params.fastq + '/*r_{1,2}.fastq.gz']

    // FastQC
    if(perform_mapping) {
        if(params.skip_fastqc){
            println("Skipping FastQC")
        } else {
            println("Running FastQC")

            Channel
            .fromPath(fastq_files)
            .set { reads_ch }

            FASTQC(reads_ch,  check_setup_ch)
        }
    }
 

    // Trimming
    if(params.trim){
        Channel
        .fromFilePairs(fastq_files)
        .set { read_pairs_ch }

        TRIMGALORE_TRIM(read_pairs_ch, check_setup_ch)
    }


    //  Concatenate FASTQ files from the same sub-library split over more than one lane
    if(params.concatenate){
        if(params.trim){
            fastq_ch = TRIMGALORE_TRIM.out.trimmed_seqs
            fastq_ch = fastq_ch.map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }.groupTuple()
        } else {
            Channel
            .fromFilePairs(fastq_files, flat: true)
            .map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }
            .groupTuple() 
            .set { fastq_ch }
        }
            CONCATENATE_FILES(fastq_ch, check_setup_ch)
    }


    //  Map the data
    if(perform_mapping) {
        println("Mapping data with split-pipe")

        if(params.concatenate) {   // Trimmed or non-trimmed concatenated FASTQ files
            read_pairs_ch = CONCATENATE_FILES.out
        } else if(params.trim){     // Trimmed but not concatenated
            read_pairs_ch = TRIMGALORE_TRIM.out
            read_pairs_ch = read_pairs_ch.map { prefix, file1, file2 -> tuple(prefix, tuple(file1, file2)) }
        } else {    // Original FASTQ files
            Channel
            .fromFilePairs(fastq_files)
            .set { read_pairs_ch }
        }
        
        if(build_index){
            index_ch = SPLITPIPE_INDEX.out
            SPLITPIPE_MAP(read_pairs_ch, index_ch, file(params.samp_list), params.chemistry, check_setup_ch)
        } else {
            SPLITPIPE_MAP(read_pairs_ch, file(params.genome_dir), file(params.samp_list), params.chemistry, check_setup_ch)
        }
    }


    // Combine sublibraries - when necessary
    //SPLITPIPE_MAP.out.count().view()
    if (perform_mapping) {
        map_ch = SPLITPIPE_MAP.out.collect().filter { v -> v.size() > 1 }   // Only combine if there is more than 1 sublibrary

        if(build_index){
            SPLITPIPE_COMBINE(map_ch, index_ch)
        } else {
            //mychannel = SPLITPIPE_MAP.out.collect()    // Check this
            SPLITPIPE_COMBINE(map_ch, file(params.genome_dir))
        }
    }


    // spipe-summary results
    if(perform_mapping) {
        println("Performing QC on split-pipe mapping")
        map_ch2 = SPLITPIPE_MAP.out.collect()   // Includes output even if only 1 sublibrary (unlike map_ch)
        SPLITPIPE_QC(map_ch2)
    }


    // Single cell QC
    if(perform_mapping) {
        splitpipe_map_ch = SPLITPIPE_MAP.out.collect()    //  Finish mapping before moving on - else failure in later step causes mapping to stop
        splitpipe_combine_ch = SPLITPIPE_COMBINE.out

        SCANPY_SINGLECELLQC(splitpipe_map_ch, "separate")
        SCANPY_SINGLECELLQC2(splitpipe_combine_ch, "combined")
    }


    // Make Seurat object
    if(perform_mapping){
        SEURAT_BUILDSEURATOBJECT(splitpipe_map_ch, "separate")
        SEURAT_BUILDSEURATOBJECT2(splitpipe_combine_ch, "combined")
    }


    // MultiQC
    if(perform_mapping){
        if(! params.skip_fastqc){
            results_ch = FASTQC.out.concat( map_ch2 )   // Include the mapping results
            results_ch = results_ch.collect()
            MULTIQC(results_ch)
        }
    }
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}


def checkParameters() {

    // Should we build the genome index file?
    if( (params.fasta == null) && (params.gtf == null) && (params.genome_name == null) && (params.genome_dir == null) ) {
        error("Error:\nSpecify a genome index folder --genome_dir, OR make a genome index denovo with --fasta/--gtf/--genome_name")
    }

    if ( (params.fasta != null) || (params.gtf != null) || (params.genome_name != null) ) {
        if ( (params.fasta == null) || (params.gtf == null) || (params.genome_name == null) ) {
            error("Error:\nParameters --fasta/--gtf/--genome_name need to be specified altogether, or not at all!")
        }
        build_index = true
    }

    //  Don't specify a new genome to build, if one already exists
    if ( (params.genome_dir != null) && (build_index == true) ) {
        error("Error:\nDon't specify an existing genome (--genome_dir) if also building a genome index from scratch!")
    }

    // Ensure an equal number of FASTA, GTF and genome name(s) have been specified
    if (build_index == true) {
       if ( (params.gtf.count(',') != params.fasta.count(',')) || (params.gtf.count(',') != params.genome_name.count(',')) ) {
            error("Error:\nParameters --fasta, --gtf and --genome_name need the same number of parameters (in a comma-separated list)")
       }
    }

    // Do we have all the relevant files for mapping?
    if ( (params.fastq != null) || (params.genome_dir != null) || (params.samp_list != null) ) {  // Mapping intended
        if ( (params.fastq == null) || (params.samp_list == null) ) {
            error("Error:\nParameters --fastq/--genome_dir/--samp_list need to be specified altogether, or not at all!")
        }
        perform_mapping = true
    }
}


def getLibraryId(filename) {
    libraryID = filename.replaceFirst(/\.s_\d$/, "")
    return libraryID
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
