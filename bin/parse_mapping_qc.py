#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Takes parse mapping file and produces graphs

#/cephfs/swingett/building_parse_pipeline/results/splitpipe/test1/process

# [sublibrary_name]/process/pre_align_stats.csv


# In[ ]:


import pandas as pd
import seaborn as sns
import glob
import os
from matplotlib import pyplot as plt


# In[ ]:


categories_to_select = pd.Series(['number_of_reads',
                                  'number_of_tscp',
                                  'reads_align_input',
                                  'reads_align_multimap',
                                  'reads_align_too_many_loci',
                                  'reads_align_unique',
                                  'reads_ambig_bc',
                                  'reads_map_transcriptome',
                                  'reads_valid_barcode'
                                 ])

outdir = 'mapping_summary_plots'
image_formats = ('png', 'svg', 'eps')

# Convert from absolute to percentage values
def convert_to_percentages(input_data):
    totals = input_data.sum(axis=1)
    input_data_pc = 100 * input_data.div(totals, axis=0)
    input_data_pc = input_data_pc.round(2)
    return(input_data_pc)
                                                         


# In[ ]:


#Import data
sublibrary_folders = glob.glob("*/")

mapping_results = pd.DataFrame()

for sublibrary_folder in sublibrary_folders:
  
    if sublibrary_folder == 'mapping_summary_plots/':   # When testing the script, don't include output file in input
        continue

    mapping_results_file = f'{sublibrary_folder}process/pre_align_stats.csv'   # Trailing '/' in sublibrary_folder
    
    if os.path.isfile(mapping_results_file):   # Check file present (ignores non-mapping results sub-folders)
        print(f'Reading in mapping data file: {mapping_results_file}')
        mapping_results_partial = pd.read_csv(mapping_results_file, index_col='statistic')
        mapping_results_partial['sublibrary'] = sublibrary_folder[:-1]   #Remove trailing /
        mapping_results = pd.concat([mapping_results, mapping_results_partial])
    
del mapping_results_partial


# In[ ]:


# Format
mapping_results = mapping_results.loc[categories_to_select, :]
mapping_results = mapping_results.pivot(columns='sublibrary', values='value')   # Wide format


# In[ ]:


# Calculate new data terms
mapping_results.loc['not_valid_barcode', :] = mapping_results.loc['number_of_reads', :] - mapping_results.loc['reads_valid_barcode', :]  - mapping_results.loc['reads_ambig_bc', :]
mapping_results.loc['unmapped', :] = mapping_results.loc['reads_align_input', :] - mapping_results.loc['reads_align_multimap', :]  - mapping_results.loc['reads_align_too_many_loci', :] - mapping_results.loc['reads_align_unique', :]
mapping_results.loc['other_genomic_location', :] = mapping_results.loc['reads_align_unique', :] - mapping_results.loc['reads_map_transcriptome', :]
mapping_results = mapping_results.rename(mapper={'number_of_tscp' : 'deduplicated_reads'})
mapping_results.loc['duplicate_copies', :] = mapping_results.loc['reads_map_transcriptome', :] - mapping_results.loc['deduplicated_reads', :]


# In[ ]:


mapping_results = mapping_results.transpose()
mapping_results = mapping_results.astype('int')


# In[ ]:


# Make output directories
for image_format in image_formats:
    outsubdir = f'{outdir}/{image_format}'
    
    if not os.path.exists(outsubdir):
        os.makedirs(outsubdir)


# In[ ]:


# Barcode plots
columns_to_select = ['reads_valid_barcode', 'reads_ambig_bc', 'not_valid_barcode']              
barcode_data = mapping_results.loc[:, columns_to_select]
barcode_data.plot(kind='bar', stacked=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

filename_base = f'valid_barcodes_count'
for image_format in image_formats:
    outfile = f'{outdir}/{image_format}/{filename_base}.{image_format}'
    print(f'Writing to {outfile}')
    plt.savefig(fname=outfile, bbox_inches='tight', pad_inches=0.5)

plt.show()


# In[ ]:


barcode_data_pc = convert_to_percentages(barcode_data)
barcode_data_pc.plot(kind='bar', stacked=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

filename_base = f'valid_barcodes_percentages'
for image_format in image_formats:
    outfile = f'{outdir}/{image_format}/{filename_base}.{image_format}'
    print(f'Writing to {outfile}')
    plt.savefig(fname=outfile, bbox_inches='tight', pad_inches=0.5)

plt.show()


# In[ ]:


# Mapping type plots
columns_to_select = ['reads_align_unique', 'reads_align_multimap', 'reads_align_too_many_loci', 'unmapped']              
align_data = mapping_results.loc[:, columns_to_select]
align_data.plot(kind='bar', stacked=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

filename_base = f'mapped_reads_count'
for image_format in image_formats:
    outfile = f'{outdir}/{image_format}/{filename_base}.{image_format}'
    print(f'Writing to {outfile}')
    plt.savefig(fname=outfile, bbox_inches='tight', pad_inches=0.5)

plt.show()


# In[ ]:


align_data_pc = convert_to_percentages(align_data)
align_data_pc.plot(kind='bar', stacked=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

filename_base = f'mapped_reads_percentages'
for image_format in image_formats:
    outfile = f'{outdir}/{image_format}/{filename_base}.{image_format}'
    print(f'Writing to {outfile}')
    plt.savefig(fname=outfile, bbox_inches='tight', pad_inches=0.5)

plt.show()


# In[ ]:


# Genomic location plots
columns_to_select = ['reads_map_transcriptome',  'other_genomic_location']              
location_data = mapping_results.loc[:, columns_to_select]
location_data.plot(kind='bar', stacked=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

filename_base = f'genomic_locations_count'
for image_format in image_formats:
    outfile = f'{outdir}/{image_format}/{filename_base}.{image_format}'
    print(f'Writing to {outfile}')
    plt.savefig(fname=outfile, bbox_inches='tight', pad_inches=0.5)

plt.show()


# In[ ]:


location_data_pc = convert_to_percentages(location_data)
location_data_pc.plot(kind='bar', stacked=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

filename_base = f'genomic_locations_percentages'
for image_format in image_formats:
    outfile = f'{outdir}/{image_format}/{filename_base}.{image_format}'
    print(f'Writing to {outfile}')
    plt.savefig(fname=outfile, bbox_inches='tight', pad_inches=0.5)

plt.show()


# In[ ]:


# Duplication plots
columns_to_select = ['deduplicated_reads',  'duplicate_copies']              
deduplication_data = mapping_results.loc[:, columns_to_select]
deduplication_data.plot(kind='bar', stacked=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

filename_base = f'deduplication_plots'
for image_format in image_formats:
    outfile = f'{outdir}/{image_format}/{filename_base}.{image_format}'
    print(f'Writing to {outfile}')
    plt.savefig(fname=outfile, bbox_inches='tight', pad_inches=0.5)

plt.show()


# In[ ]:


deduplication_data_pc = convert_to_percentages(deduplication_data)
deduplication_data_pc.plot(kind='bar', stacked=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

filename_base = f'deduplication_plots_pc'
for image_format in image_formats:
    outfile = f'{outdir}/{image_format}/{filename_base}.{image_format}'
    print(f'Writing to {outfile}')
    plt.savefig(fname=outfile, bbox_inches='tight', pad_inches=0.5)

plt.show()


# In[ ]:


print('Done')

