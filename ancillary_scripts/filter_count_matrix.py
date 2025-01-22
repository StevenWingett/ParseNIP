# Setup options
import argparse
import os
#import pprint

parser = argparse.ArgumentParser()
parser.add_argument("--samples", action="store", type=str, metavar='',
                    help="Path to file listing samples to extract [default=samplesheet.txt]", 
                    default='samplesheet.txt')
parser.add_argument("--metadata", action="store", type=str, metavar='',
                    help="Path to file listing the cell metadata to extract [default=cell_metadata.csv]", 
                    default='cell_metadata.csv')
parser.add_argument("--counts", action="store", type=str, metavar='',
                    help="Path count matrix file [default=count_matrix.mtx]", 
                    default='count_matrix.mtx')
parser.add_argument("--outdir", action="store", type=str, metavar='',
                    help="Path to output directory [default=subset_data]", 
                    default='subset_data')


# Parse arguments and setup an output directory
args = parser.parse_known_args()    #Use parse_known_arg to differentiate between arguments pre-specified and those that are not
options = args[0]   # Get the 2 arrays of known/unknown arguments from the tuple
#files = set(args[1])


# Make the output directory
if os.path.exists(options.outdir):
    sys.exit(f'Output directory {"options.outdir"} already exists, adjust --outdir option')
    pass
else:
    os.makedirs(options.outdir)


# Identify samples of interest
samples_of_interest = {}

print(f'Reading samples of interest file "{options.samples}"')

i = 0
with open(options.samples, 'r') as f: 
    for line in f:
        i += 1
        samples_of_interest[line.rstrip()] = None    # .rstrip() to remove newline characters

    
f.close()    # Close the file

print(f'{i} line(s) read and {len(samples_of_interest)} unique entry/entries found')
#pprint.pprint(samples_of_interest)


# Read cell metadata and extract cells derived from the samples of interest
# Create an object that stores cell metadata old line numbers (minus 1) and new line numbers
# Record genes (as line numbers of interest)
print(f'Reading cell metadata file "{options.metadata}"')
i = 0
valid_cells_count = 0
cell_metadata_lookup = {}    # {position_in_old_file} = position_in_new_file

with open(options.metadata, 'r') as f, open(f'{options.outdir}/cell_metadata.csv', 'w') as fout:

    fout.write(f.readline())   # Write out header

    for line in f:
        i += 1
        sample = line.split(',')[1]

        if sample in samples_of_interest:
            valid_cells_count += 1
            cell_metadata_lookup[i] = valid_cells_count
            fout.write(line)
             
# Close the files
f.close()
fout.close()

print(f'{valid_cells_count} cells extracted out of a total of {i}')


# Read the count matrix file
print(f'Reading the count matrix file "{options.counts}"')
i = 0
valid_cell_gene_count_observations = 0
#total_transcript_count = 0
#cell_metadata_lookup = {}    # {position_in_old_file} = position_in_new_file

with open(options.counts, 'r') as f, open(f'{options.outdir}/_tmp_count_matrix.mtx', 'w') as fout:

    f.readline()   # Skip headers
    f.readline()
    gene_count = f.readline().split(' ')[1]   # Use later when making header

    for line in f:
        i += 1
        cell_id, gene_id, transcript_count = line.rstrip().split(' ')
        cell_id = int(cell_id)
        gene_id = int(gene_id)
        transcript_count = int(transcript_count)

        #print(cell_id)

        if cell_id in cell_metadata_lookup:
            valid_cell_gene_count_observations += 1
            new_cell_id = cell_metadata_lookup[cell_id]
            edited_line = f'{new_cell_id} {gene_id} {transcript_count}\n'
            fout.write(edited_line)
            #total_transcript_count += transcript_count
             
# Close the files
f.close()
fout.close()

print(f'{valid_cell_gene_count_observations} valid cell/gene count observations extracted out of a total of {i}')


# Make the final output file
# First make header
header = '%%MatrixMarket matrix coordinate integer general\n'
header += f'%Rows=cells ({valid_cells_count}), Cols=genes ({gene_count})\n'
header += f'{valid_cells_count} {gene_count} {valid_cell_gene_count_observations}\n'


with open(f'{options.outdir}/_header', 'w') as fout:
    fout.write(header)
fout.close()  # Close the file

# Then append the body of the file
command = f'cat {options.outdir}/_header {options.outdir}/_tmp_count_matrix.mtx > {options.outdir}/count_matrix.mtx'
os.system(command)
os.remove(f'{options.outdir}/_header')
os.remove(f'{options.outdir}/_tmp_count_matrix.mtx')

print(f'Writing output to "{options.outdir}"')

print('Done')
