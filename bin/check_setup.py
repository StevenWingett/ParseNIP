import os
import glob
import pprint
import re


# Get arguments
#fastq_folder = "FASTQ"
fastq_folder = "FASTQTEST"
genome_dir = 'GENOME_INDEX'
fasta_file = 'a.txt'
gtf_file = 'b.txt'
samplesheet_file = 'parse_samplesheet.txt'


# # # # # # # # # # # # # # # # # # # #
# Functions
# # # # # # # # # # # # # # # # # # # #

def check_fastq(fastq_folder):   # Is the FASTQ folder ok?

    problem_flag = 0

    # Is FASTQ folder present?
    if not os.path.isdir(fastq_folder):
        print(f'{fastq_folder} not found!')
        problem_flag = 1

    # Are FASTQ files present?
    fastq_files = glob.glob(f'{fastq_folder}/*.fastq.gz')
    fastq_files = fastq_files + glob.glob(f'{fastq_folder}/*.fq.gz')              

    if len(fastq_files) == 0:
        print(f'No FASTQ files found in {fastq_folder}')
        problem_flag = 1
    
    # Do all the FASTQ files have the correct file structure?
    # [identifier].s_[lane number].r_[read1 or read2].fq.gz
    for i in range(len(fastq_files)):

        # Let's simplify the filenames first
        fastq_files[i] = os.path.basename(fastq_files[i])

        if fastq_files[i][-6:] == '.fq.gz':
            fastq_files[i] = fastq_files[i][:-6]
        else:
            fastq_files[i] = fastq_files[i][:-9]

        pat = '.+\.s_\d+\.r_[12]$'

        pat_match = re.search(pat, fastq_files[i])

        if not pat_match:
            print(f'Incorrect FASTQ filename structure : {fastq_files[i]}')
            problem_flag = 1

    # Create lookup dictionary
    fastq_files_dict = {}

    for fastq_file in fastq_files:   
        fastq_files_dict[fastq_file] = None

    # *.fastq.gz same as *.fq.gz e.g. test1.s_1.r_1.fastq.gz AND test1.s_1.r_1.fq.gz
    if len(fastq_files) != len(fastq_files_dict):
        print('Some files have the same [identifier].s_[lane number].r_[read1 or read2] structure!')
        print('Is there a *.fastq.gz with the same name as a *.fq.gz file?')
        print('E.g.test1.s_1.r_1.fastq.gz AND test1.s_1.r_1.fq.gz')
        problem_flag = 1

    # Is each FASTQ file a member of a pair?
    for fastq_file in fastq_files:   # Now check
        if fastq_file[-1] == '1':
            lookup_file = fastq_file[:-1] + '2'
        else:
            lookup_file = fastq_file[:-1] + '1'

        if lookup_file not in fastq_files_dict:
            print(f'{lookup_file} has no matching pair!')
            problem_flag = 1
  
    return problem_flag



def check_genome_dir(genome_dir):    # Is the pre-build genome ok?

    problem_flag = 0

    # Is genome directory folder present?
    if not os.path.isdir(genome_dir):
        print(f'{genome_dir} not found!')
        problem_flag = 1
        return problem_flag

    # Check if one of the expected files is present
    expected_file = 'all_genes.csv'
    expected_file = f'{genome_dir}/{expected_file}'

    if not os.path.isfile(expected_file):
        print(f'{expected_file} not found!')
        problem_flag = 1

    return problem_flag



# check_samplesheet
# Samplesheet should be space-delimited file, listing each sample and well ID e.g.
# Sample1 A1
# Sample2 A2
# .
# .
def check_samplesheet(samplesheet_file):
    problem_flag = 0
    allowed_wells_dict = {}   # Valid plate well IDs

    for char in 'ABCDEFGH':
        for i in range(1, 13):
            well_id = char + str(i)
            allowed_wells_dict[well_id] = 0    # Initialise at 0

    with open(samplesheet_file, 'r') as f:

        for line in f:
            line = line.rstrip()

            # Check not only whitespace characters
            if line == '':
                print(f'Empty line in {samplesheet_file}!')
                return(1)
            
            pat_match = re.search('^\s*$', line)
            if(pat_match):
                print(f'Empty line in {samplesheet_file}!')
                return(1)

            # Not exactly two elements after splitting on whitespace
            line_elements = line.split(' ')
            if len(line_elements) != 2:
                print(f'Line does not contain exactly 2 elements:\n{line}')
                return(1)

            # Not a valid well ID
            if line_elements[1] not in allowed_wells_dict:
                print(f'Invalid Well ID: {line_elements[1]}')
                return(1)

            # Samples should only be alpha-numeric characters
            pat_match = re.search('^[A-z0-9_]+$', line_elements[0])
            if not pat_match:
                print(f'Sample Names need to only be alpha-numeric characters, NOT {line_elements[0]}!')
                return(1)

            # Well used more than once
            allowed_wells_dict[line_elements[1]] = allowed_wells_dict[line_elements[1]] + 1
            if allowed_wells_dict[line_elements[1]] > 1:
                print(f'Well ID {allowed_wells_dict[line_elements[1]]} listed more than once!')
                return 1
            
    return 0
            



# # # # # # # # # # # # # # # # # # # #
# Main Code
# # # # # # # # # # # # # # # # # # # #

#  Check FASTQ folder
exit_code = check_fastq(fastq_folder)


# Check FASTA file, if specified
if not os.path.isfile(fasta_file):
    print(f'FASTA file {fasta_file} not found!')
    exit_code += 1


# Check GTF file, if specified
if not os.path.isfile(gtf_file):
    print(f'GTF file {gtf_file} not found!')
    exit_code += 1

# Check Genome folder, if specified
exit_code = exit_code + check_genome_dir(genome_dir)


# Check Samplesheet
print(samplesheet_file)
exit_code = exit_code + check_samplesheet(samplesheet_file)


if exit_code == 0:
    print('Done')
else: 
    print('Problem(s) with setup!')

exit(exit_code)


