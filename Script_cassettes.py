"""
Python Script
To: Build cassettes for each genome_id__contig_id
Inputs: 
'test6_allgenomes_out-DGR_example.gtf'
'test6-09_nrseqRT_names_example'  #File with the names of the 1012 nrRT sequences
"""

import pandas as pd
import numpy as np

# 1. Read test6_allgenomes_out-DGR.gtf
path_fullgtf = './test6_allgenomes_out-DGR.gtf'                     #1 ARGS   #the _example.gtf has 5 genome_ids
fullgtf = pd.read_csv(path_fullgtf, sep = '\t') 
#print(fullgtf[['genome_id', 'count_DGRcomp', '#ID', 'DGRcomp', 'string', 'original_start', 'original_end', 'start', 'end', 'A-to-N-subs', 'non-A-to-N-subs']])
#print(fullgtf[['genome_id', 'count_DGRcomp', '#ID', 'DGRcomp', 'string', 'original_start', 'original_end', 'start', 'end']])

# 2. Replace asterisks '*' with 'start' and 'end' columns
def replace_asterisk(row):
    if row['original_start'] == '*':
        row['original_start'] = row['start']
    if row['original_end'] == '*':
        row['original_end'] = row['end']
    return row        

fullgtf = fullgtf.apply(lambda row: replace_asterisk(row), axis = 1) #axis = 1 indica fila
#print(fullgtf[['genome_id', 'count_DGRcomp', '#ID', 'DGRcomp', 'string', 'original_start', 'original_end', 'start', 'end']])

# 3. For each genome_id__contig_id, keep the rows with 'TR' that have unique 'original_start' and 'original_end' (i.e. drop duplicates)
# 4. Reorder TR, VR, RT based on column 'original_start'
# 5. Count RTs 
# 6. Count TRs
# 7. Create columns 'check_cassette' and 'order_cassette'

new_fullgtf = []
df_cassette = {'genome_id':[], 'contig_id': [],'check_cassette': [], 'order_cassette':[]}

for keys, df_keys in fullgtf.groupby(['genome_id', '#ID']):
    #print(keys, '\n', df_keys, '\n')
    df_keys_TR = df_keys[df_keys['DGRcomp'] == 'TR']            #the df_keys_TR dataframe has the rows with TR
    df_keys = df_keys[df_keys['DGRcomp'] != 'TR']               #the df_keys dataframe has the rows with VR and RT
    df_keys_TR = df_keys_TR.drop_duplicates(['original_start', 'original_end'])
    #print(df_keys_TR['original_start'].dtypes)  # <-- output: object
    #print(df_keys_TR['original_end'].dtypes)  # <-- output: object
    #df_keys_TR['original_start'] = df_keys_TR['original_start'].astype(float) # <--
    #df_keys_TR['original_end'] = df_keys_TR['original_end'].astype(float) # <--
    #df_keys_TR = df_keys_TR.astype({'original_start': int}) # <--
    #df_keys_TR = df_keys_TR.astype({'original_end': int}) # <--
    
    df_keys = pd.concat([df_keys_TR, df_keys])
    df_keys = df_keys.astype({'original_start': int}) # <--
    df_keys = df_keys.astype({'original_end': int}) # <--
    df_keys = df_keys.sort_values(['original_start', 'original_end']) #ascending = True by default 
    #df_keys = df_keys.sort_values('original_start') #ascending = True by default   # <--

    count_RT = len(df_keys[df_keys['DGRcomp'] == 'RT'])
    df_keys['count_RT_percontig'] = count_RT
    count_TR = len(df_keys[df_keys['DGRcomp'] == 'TR'])

    check_cassette = None
    if count_RT >= 2:
        check_cassette = 'CHECK_2RT'
    if count_TR >= 2:
        check_cassette = 'CHECK'
    
    order_cassette = []
    for dgrcomp, symbol in zip(df_keys['DGRcomp'], df_keys['string']):
        order_cassette.append(symbol+dgrcomp)       #example: +VR
    order_cassette = ','.join(order_cassette)       #example: +VR,+TR,+RT

    genome_id = df_keys['genome_id'].values[0]
    contig_id = df_keys['#ID'].values[0]

    df_cassette['genome_id'].append(genome_id)
    df_cassette['contig_id'].append(contig_id)
    df_cassette['check_cassette'].append(check_cassette)
    df_cassette['order_cassette'].append(order_cassette)
    
    new_fullgtf.append(df_keys)

new_fullgtf = pd.concat(new_fullgtf)
#print(new_fullgtf[['genome_id', 'count_DGRcomp', 'count_RT_percontig', '#ID', 'DGRcomp', 'string', 'original_start', 'original_end', 'start', 'end']])
#new_fullgtf.to_csv('./test6_allgenomes_out-DGR_ordered.gtf', sep = '\t')           #1 ARGS #This is a new file 'allgenomes_out-DGR.gtf', ordered by column 'original_start', and there is a new column 'count_RT_percontig'

df_cassette = pd.DataFrame(df_cassette)
#print(df_cassette)
#df_cassette.to_csv('./test6_order_cassette_2263.csv' , sep = ';')                  #2 ARGS  #genome_id, contig_id, check_cassette, order_cassette, genome__contig  #It should have 2232 rows  (I think they correspond to the unique contig_ids)
#df_cassette.to_csv('./final_outputs/test6_order_cassette_2263.csv' , sep = ';')    #2 ARGS

# 8. Show unique values in column 'order_cassette' + count
df_count = df_cassette[['order_cassette']].groupby('order_cassette').size()
#print(df_count)
#df_count.to_csv('./test6_order_cassette_2263_count.csv', sep = ';')                #3 ARGS
#df_count.to_csv('./final_outputs/test6_order_cassette_2263_count.csv', sep = ';')  #3 ARGS

#=======================================================================
# Now I want to know which are the 1012 non-redundant genome_ids__num_RT.
# But there is a problem: 
#   The file './test6_order_cassette_2263.csv' has the 'genome_id__contig_id'
#   And the file 'test6-09_nrseqRT' has the headers as 'genome_id__num_RT'

# In Sophon:
# less test6-09_nrseqRT | grep '>' > test6-09_nrseqRT_names  #>3300025682_87__1

# 9. Read the list of genome_ids__contig_ids from test6-09_nrseqRT_names
path_listnames = './test6-09_nrseqRT_names'                                     #4 ARGS  #3300025682_87__1
listnames = pd.read_csv(path_listnames, header = None)

listnames = listnames.rename({0:'genome__numRT'}, axis = 1)
listnames['genome_id'] = listnames['genome__numRT'].str.split('__').str[0]  #genome_id__num_RT
listnames['num_RT'] = listnames['genome__numRT'].str.split('__').str[1]     #genome_id__num_RT
#print(listnames) #3columns: genome__numRT, genome_id, num_RT  #It should have 1012 nrRT rows

# 10. Back to the file 'new_fullgtf', add the 'num_RT' per genome_id
new_fullgtf_RT = new_fullgtf[new_fullgtf['DGRcomp'] == 'RT']
new_fullgtf_RT = new_fullgtf_RT.reset_index()
#print(new_fullgtf_RT[['genome_id', '#ID', 'DGRcomp', 'count_RT_percontig']])

new_fullgtf_RT_num = []
for keys, df_genomeid in new_fullgtf_RT.groupby('genome_id'):
    #print(keys, '\n', df_genomeid, '\n')
    df_genomeid['num_RT'] = np.arange(1, len(df_genomeid) +1)
    new_fullgtf_RT_num.append(df_genomeid)      #append the dataframes to my list of dataframes

new_fullgtf_RT_num = pd.concat(new_fullgtf_RT_num)     #concat the list of dataframes
new_fullgtf_RT_num['genome__numRT'] = new_fullgtf_RT_num[['genome_id', 'num_RT']].astype(str).agg('__'.join, axis = 1) #https://www.statology.org/pandas-combine-two-columns/
#print(new_fullgtf_RT_num[['genome_id', '#ID', 'DGRcomp', 'count_RT_percontig', 'num_RT', 'genome__numRT']]) #It should have 2263 RT rows

# 11. Merge two dataframes 'filenames' (which has the 1012 nrRT sequence names) and 'new_fullgtf_RT_num' 
new_fullgtf_RT_num = new_fullgtf_RT_num[['#ID', 'genome__numRT']]
listnames = listnames.merge(new_fullgtf_RT_num, on = 'genome__numRT', how = 'left')

listnames['genome__contig'] = listnames[['genome_id', '#ID']].astype(str).agg('__'.join, axis = 1)
#print(listnames) #It should have 1012 nrRT rows

#===============================================================================
# 12. Back to the dataframe 'df_cassette', I will add a new column 'genome_contig'
df_cassette['genome__contig'] = df_cassette[['genome_id', 'contig_id']].astype(str).agg('__'.join, axis = 1)
#print(df_cassette)    #It should have 2232 rows  (I think they correspond to the unique contig_ids)

# 13. Keep only rows that are on the list
list_1012 = listnames['genome__contig'].to_list()
df_cassette_1012 = df_cassette[df_cassette['genome__contig'].isin(list_1012)]
#print(df_cassette_1012)
#df_cassette_1012.to_csv('./test6_order_cassette_1012.csv', sep = ';')                          #5 ARGS
#df_cassette_1012.to_csv('./final_outputs/test6_order_cassette_1012.csv', sep = ';')            #5 ARGS  
                                         
# 14. Show unique values in column 'order_cassette' + count
df_count_1012 = df_cassette_1012[['order_cassette']].groupby('order_cassette').size().sort_values(ascending = False)
#print(df_count_1012)
#df_count_1012.to_csv('./test6_order_cassette_1012_count.csv', sep = ';')                       #6 ARGS
#df_count_1012.to_csv('./final_outputs/test6_order_cassette_1012_count.csv', sep = ';')         #6 ARGS

#=============================================================================================
# 15. I manually checked the file 'test6_order_cassette_1012.csv' and saved it as manual. 
# Then, show unique values in column 'order_cassette' + count
#df_cassette_1012_manual = pd.read_csv('./test6_order_cassette_1012_manual.csv', sep = ';')                 #7 ARGS
df_cassette_1012_manual = pd.read_csv('./final_outputs/test6_order_cassette_1012_manual.csv', sep = ';')   #7 ARGS

df_count_1012_manual = df_cassette_1012_manual[['order_cassette']].groupby('order_cassette').size().sort_values(ascending = False)
#df_count_1012_manual.to_csv('./test6_order_cassette_1012_manual_count.csv', sep = ';')                     #8 ARGS
df_count_1012_manual.to_csv('./final_outputs/test6_order_cassette_1012_manual_count.csv', sep = ';')       #8 ARGS
