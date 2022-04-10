"""
Python script
Inputs:
- 'test6_funannotation_complete_VRposition.csv'
- input2 -->  directory './aaseq' with all the 'product.tsv' files (output from InterproScan)
Output:
- new directory './aaseq_complete' with all the 'product_complete.tsv' files (merged columns, reset coordinates, and column 'overlap_domains')
"""

import pandas as pd
from os import listdir, makedirs
from os.path import isfile, join, exists
import numpy as np

path_funannot = './test6_funannotation_complete_VRposition.csv'    #1 INPUT
funannot = pd.read_csv(path_funannot, sep = ';')
funannot = funannot[['ORF_id', 'ORF_start', 'ORF_end', 'VR_start', 'VR_end', 'VR_position', 'VR_strand', 'biome', 'ecosystem_category', 'product']]

count = funannot['ORF_id'].value_counts().to_dict()
product = funannot.set_index('ORF_id').to_dict()            # <-- 
#print(product)
dict_comb = {key:[count[key], product[key]] for key in count}  # <--
print(dict_comb)

for orfid in count:
    value = count[orfid]
    if value > 1:
        print(f'The ORF_id {orfid} is repeated, y the product is')  # <--

#*For each file 'product.tsv' in the directory 'aaseq'
#*Add column names
#*Search each unique 'ORF_id' on the file 'test6_funannotation_complete_VRposition.csv'
#*If 'ORF_id' was only found once, merge by the column 'ORF_id' with 9 columns: 'ORF_start', 'ORF_end', 'VR_start', 'VR_end', 'VR_position', 'VR_strand', 'biome', 'ecosystem_category', 'product'
#*If 'ORF_id' was found multiple times,  print message with the name of the 'ORF_id', and duplicate domain annotations

path_aaseq = './aaseq/'                                                                     #2 INPUT
files_aaseq = [join(path_aaseq, f) for f in listdir(path_aaseq) if isfile(join(path_aaseq, f)) and f[-3:] == 'tsv'] #and '_complete' not in f

if not exists('aaseq_complete'):
    makedirs('aaseq_complete')

# Add column names:       
columnsnames = ['ORF_id', 'seq_md5', 'seq_length', 'analysis', 'signature_id', 'signature_desc', 'dom_start', 'dom_end', 'score', 'status', 'date', 'interpro_id', 'interpro_desc']
for path in files_aaseq:
    df = pd.read_csv(path, sep = '\t', header = None)
    df.columns = columnsnames
    # Merge by the column 'ORF_id':
    df = pd.merge(df, funannot, on = 'ORF_id', how = 'left')
    #df.to_csv('./archivode_verificacion.csv', sep = ';')
    #count = df['ORF_id'].value_counts().to_dict()       ### I shouldn't count here
    #for orfid in count:         
        #value = count[orfid]
        #if value > 1: 
            #print(f'The ID {orfid} is repeated in the file {path}')
    originalpath = path.split('/')  #['.', 'aaseq', 'atp_1.tsv']
    newname = originalpath[-1].split('.')   #['atp_1_complete', 'tsv']
    newname[0] += '_complete'   #newname[0] + '_complete'
    newname = '.'.join(newname) #atp_1_complete.tsv
    originalpath[-1] = newname  #['.', 'aaseq', 'atp_1_complete.tsv']
    originalpath[-2] += '_complete'
    originalpath = '/'.join(originalpath) #./aaseq/atp_1_complete.tsv
    #print(originalpath)
    #df.to_csv(originalpath, index = False, sep = ';')  #originalpath is a string

#*Then, reset the numbers in 'ORF_start', 'ORF_end', 'VR_start', 'VR_end' to their equivalents starting from 1 (('ORF_start' - 'ORF_start +1) and ('ORF_end' 'VR_start' 'VR_end' - 'ORF_start' + 1)/3)
#*Put those values into new columns 'ORF_newstart', 'ORF_newend', 'VR_newstart', 'VR_newend'
#*Don't forget to keep the values as int and not as float.
#*Save these tsv files (just in case)
#*Drop columns 'ORF_start', 'ORF_end', 'VR_start', 'VR_end'
#*Save these tsv files

    # Then, reset the numbers in 'ORF_start', 'ORF_end', 'VR_start', 'VR_end':
    df['ORF_newstart'] = df['ORF_start'] - df['ORF_start'] + 1
    df['ORF_newend'] = np.ceil((df['ORF_end'] - df['ORF_start'] + 1)/3)
    df['VR_newstart'] = np.ceil((df['VR_start'] - df['ORF_start'] + 1)/3)
    df['VR_newend'] = np.ceil((df['VR_end'] - df['ORF_start'] + 1)/3)
    #df.to_csv(originalpath, index = False, sep = ';')

    #df.drop(columns=['B', 'C']) o df.drop(['B', 'C'], axis=1)
    df.drop(columns=['ORF_start', 'ORF_end', 'VR_start', 'VR_end']) 
    #df.to_csv(originalpath, index = False, sep = ';')

#*Next, create column 'overlap_domain'
#*If 'VR_newstart' is within 'start' and 'end' or if 'VR_newend' is within 'start' and 'end', put YES on 'overlap_domain'
#*Else, put NO
#*Filtrate rows with YES    

    # Next, create column 'overlap_domain':
    df.loc[(df['dom_start'] < df['VR_newstart']) & (df['VR_newstart'] < df['dom_end']), 'overlap_domain'] = 'YES'
    df.loc[(df['dom_start'] < df['VR_newend']) & (df['VR_newend'] < df['dom_end']), 'overlap_domain'] = 'YES'
    df['overlap_domain'] = df['overlap_domain'].fillna('NO')

    # Save file:
    # 13 columns in .tsv & 9+4-4 columns in funannotation.csv (don't consider 'seq_md5', 'date', 'ORF_newstart', 'ORF_newend')
    final_columnsnames = ['ORF_id', 'seq_length', 'analysis', 'signature_id', 'signature_desc', 'dom_start', 'dom_end', 'VR_newstart', 'VR_newend', 'overlap_domain', 'VR_position', 'VR_strand', 'score', 'status', 'interpro_id', 'interpro_desc', 'biome', 'ecosystem_category', 'product']
    df = df[final_columnsnames]
    df.to_csv(originalpath, index = False, sep = ';')           #3 OUTPUT
