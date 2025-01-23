#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta = getGenomeAttribute('fasta')
params.gtf = getGenomeAttribute('gtf')
params.genome_name = null
params.genome_dir = null
params.chemistry = 'v3'
params.samp_list = null
params.concatenate = false
params.fastq = null
params.trim = false
params.genome = null

log.info """\
    P A R S E  E V E R C O D E - N F  P I P E L I N E
    =================================================
    Output folder                   : ${params.outdir}
    FASTA folder                    : ${params.fasta}
    GTF file                        : ${params.gtf}
    Genome Name                     : ${params.genome_name}
    Genome                          : ${params.genome}
    Pre-built genome index folder   : ${params.genome_dir}
    FASTQ files folder              : ${params.fastq}
    Chemistry version               : ${params.chemistry}
    Sample list file                : ${params.samp_list}
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

include { SPLITPIPE_INDEX } from './modules/local/splitpipe/index/main.nf'
include { SPLITPIPE_MAP } from './modules/local/splitpipe/map/main.nf'
include { CONCATENATE_FILES } from './modules/local/cat/contatenate_files/main.nf'
include { SPLITPIPE_COMBINE } from './modules/local/splitpipe/combine/main.nf'
include { TRIMGALORE_TRIM } from './modules/local/trim_galore/trim_galore_hardtrim/main.nf'
include { FASTQC } from './modules/local/fastqc/fastqc/main.nf'
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

     def fastq_files = [params.fastq + '/*r_{1,2}.fq.gz', params.fastq + '/*r_{1,2}.fastq.gz']

    // FastQC
    Channel
    .fromPath(fastq_files)
    .set { reads_ch }

    FASTQC(reads_ch)
 

    // Trimming
    if(params.trim){
        Channel
        .fromFilePairs(fastq_files)
        .set { read_pairs_ch }

        TRIMGALORE_TRIM(read_pairs_ch)
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
            CONCATENATE_FILES(fastq_ch)
    }


    // Create split-pipe index file
    if(build_index) {
        println("Building genome index")
        SPLITPIPE_INDEX(file(params.fasta), file(params.gtf), params.genome_name)
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
            SPLITPIPE_MAP(read_pairs_ch, index_ch, file(params.samp_list), params.chemistry)
        } else {
            SPLITPIPE_MAP(read_pairs_ch, file(params.genome_dir), file(params.samp_list), params.chemistry)
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


    // Mapping QC
    // TODO

    // spipe-summary results
    // TODO

    // Make Seuarat / Scanpy objects
    // TODO

    // Single cell QC
    // TODO


    // MultiQC
    results_ch = FASTQC.out.collect()
    MULTIQC(results_ch)
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
    if ( (params.genome_dir != null) && (build_index == true) ){
        error("Error:\nDon't specify an existing genome (--genome_dir) if also building a genome index from scratch!")
    }

    // Do we have all the relevant files for mapping?
    if ( (params.fastq != null) || (params.genome_dir != null) || (params.samp_list != null) ) {  // Mapping intended
        if ( (params.fastq == null) || (params.samp_list == null) ) {
            error("Error:\nParameters --fastq/--samp_list need to be specified altogether, or not at all!")
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
