"""
Python Script 
To: Build cassettes for each genome_id__contig_id
Inputs & Outputs: 
input1 --> 'test6_allgenomes_out-DGR_example.gtf'                   #1
output1 --> './test6_order_cassette_2263_args.csv'                  #2
output2 --> './test6_order_cassette_2263_count_args.csv'            #3
input2 --> 'test6-09_nrseqRT_names_example'                         #4  #File with the names of the 1012 nrRT sequences
output3 --> './test6_order_cassette_1012.csv'                       #5
output4 --> './test6_order_cassette_1012_count.csv'                 #6  

# Example of how to run the script:
#python .\script_cassettes.py -input1 'test6_allgenomes_out-DGR_example.gtf' -output1 'out1.csv' -output2 'out2.csv' -input2 'test6-09_nrseqRT_names_example2' -output3 'out3.csv' -output4 'out4.csv'
"""

import argparse             #Step1
import pandas as pd
import numpy as np

def main(args):
    path_fullgtf = args.input1
    path_output1 = args.output1
    path_output2 = args.output2
    path_listnames = args.input2
    path_output3 = args.output3
    path_output4 = args.output4


    # 1. Read test6_allgenomes_out-DGR.gtf
    fullgtf = pd.read_csv(path_fullgtf, sep = '\t')                      #1 ARGS   #path_fullgtf = './test6_allgenomes_out-DGR.gtf'     #the _example.gtf has 5 genome_ids
    #print(fullgtf[['genome_id', 'count_DGRcomp', '#ID', 'DGRcomp', 'string', 'original_start', 'original_end', 'start', 'end', 'A-to-N-subs', 'non-A-to-N-subs']])
    #print(fullgtf[['genome_id', 'count_DGRcomp', '#ID', 'DGRcomp', 'string', 'original_start', 'original_end', 'start', 'end']])

    # 2. Replace asterisks '*' with 'start' and 'end' columns
    def replace_asterisk(row):
        if row['original_start'] == '*':
            row['original_start'] = row['start']
        if row['original_end'] == '*':
            row['original_end'] = row['end']
        return row        

    fullgtf = fullgtf.apply(lambda row: replace_asterisk(row), axis = 1) #axis = 1 indicates row
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
        
        df_keys = pd.concat([df_keys_TR, df_keys])
        df_keys = df_keys.astype({'original_start': int}) # <-
        df_keys = df_keys.astype({'original_end': int}) # <-
        df_keys = df_keys.sort_values(['original_start', 'original_end']) #ascending = True by default 
        #df_keys = df_keys.sort_values('original_start') #ascending = True by default   # <-

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
    #new_fullgtf.to_csv('./test6_allgenomes_out-DGR_ordered.gtf', sep = '\t')                #1extra ARGS #This is a new file 'allgenomes_out-DGR.gtf', ordered by column 'original_start', and there is a new column 'count_RT_percontig'

    df_cassette = pd.DataFrame(df_cassette)
    df_cassette.to_csv(path_output1 , sep = ';')                                             #2 ARGS  #df_cassette.to_csv('./test6_order_cassette_2263_args.csv' , sep = ';')   #It should have 2232 rows  (I think they correspond to the unique contig_ids)

    # 8. Show unique values in column 'order_cassette' + count
    df_count = df_cassette[['order_cassette']].groupby('order_cassette').size()
    df_count.to_csv(path_output2, sep = ';')                                                 #3 ARGS  #df_count.to_csv('./test6_order_cassette_2263_count_args.csv', sep = ';')

#==========================================================================================================================================================================================
    # Now I want to know which are the 1012 non-redundant genome_ids__numRT.
    # But there is a problem: 
       #The file './test6_order_cassette_2263.csv' has the 'genome_id__contig_id'
       #And the file 'test6-09_nrseqRT' has the headers as 'genome_id__numRT'

    # In Sophon:
    # less test6-09_nrseqRT | grep '>' > test6-09_nrseqRT_names  #>3300025682_87__1       

    # 9. Back to the dataframe 'new_fullgtf', for each unique 'genome_id', count the numRT if the row is 'RT'  
    new_fullgtf['numRT'] = 0       #or ''    #add an empty column 'numRT' per genome_id
    new_fullgtf_numRT = []
    for keys, df_genomeid in new_fullgtf.groupby('genome_id'):
        #print(keys, '\n', df_genomeid[['genome_id', '#ID', 'DGRcomp']], '\n')
        df_genomeid.loc[df_genomeid['DGRcomp'] == 'RT', ['numRT']] = np.arange(1, len(df_genomeid[df_genomeid['DGRcomp'] == 'RT']) +1)
        df_genomeid = df_genomeid.astype({'numRT': int})

        new_fullgtf_numRT.append(df_genomeid)      #append the dataframes to my list of dataframes
    
    new_fullgtf_numRT = pd.concat(new_fullgtf_numRT)     #concat the list of dataframes 
    new_fullgtf_numRT = new_fullgtf_numRT.reset_index()  
    #print(new_fullgtf_numRT[['genome_id', '#ID', 'DGRcomp', 'numRT']])                                 #1extra ARGS    #new_fullgtf_numRT[['genome_id', '#ID', 'DGRcomp', 'numRT']].to_csv(path_output1extra, sep = '\t')

    # Add a new column 'genome__numRT', combining two columns
    new_fullgtf_numRT['genome__numRT'] = new_fullgtf_numRT[['genome_id', 'numRT']].astype(str).agg('__'.join, axis = 1) 
    #print(new_fullgtf_numRT[['genome_id', '#ID', 'DGRcomp', 'count_RT_percontig', 'numRT', 'genome__numRT']]) #It should have 2263 RT rows

    # 10. Read the list of genome_ids__contig_ids from test6-09_nrseqRT_names
    if args.input2 != 'None':    #if the file exists, run this part of the script
        listnames = pd.read_csv(path_listnames, header = None)                                          #4 ARGS  #3300025682_87__1   #path_listnames = './test6-09_nrseqRT_names'

        listnames = listnames.rename({0:'genome__numRT'}, axis = 1)
        listnames['genome_id'] = listnames['genome__numRT'].str.split('__').str[0]  #genome_id__numRT
        listnames['numRT'] = listnames['genome__numRT'].str.split('__').str[1]     #genome_id__numRT
        #print(listnames) #3columns: genome__numRT, genome_id, numRT  #It should have 1012 nrRT rows

        # 11. Merge two dataframes 'filenames' (which has the 1012 nrRT sequence names) and 'new_fullgtf_numRT' 
        new_fullgtf_numRT = new_fullgtf_numRT[['#ID', 'genome__numRT']]
        listnames = listnames.merge(new_fullgtf_numRT, on = 'genome__numRT', how = 'left')

        listnames['genome__contig'] = listnames[['genome_id', '#ID']].astype(str).agg('__'.join, axis = 1)
        #print(listnames) #It should have 1012 nrRT rows

        # 12. Back to the dataframe 'df_cassette', I will add a new column 'genome_contig'
        df_cassette['genome__contig'] = df_cassette[['genome_id', 'contig_id']].astype(str).agg('__'.join, axis = 1)
        #print(df_cassette)    #It should have 2232 rows  (I think they correspond to the unique contig_ids)

        # 13. Keep only rows that are on the list
        list_1012 = listnames['genome__contig'].to_list()
        df_cassette_1012 = df_cassette[df_cassette['genome__contig'].isin(list_1012)]
        df_cassette_1012.to_csv(path_output3 , sep = ';')                                               #5 ARGS     #df_cassette_1012.to_csv('./final_outputs/test6_order_cassette_1012.csv', sep = ';')

        # 14. Show unique values in column 'order_cassette' + count
        df_count_1012 = df_cassette_1012[['order_cassette']].groupby('order_cassette').size().sort_values(ascending = False)
        df_count_1012.to_csv(path_output4 , sep = ';')                                                  #6 ARGS      #df_count_1012.to_csv('./final_outputs/test6_order_cassette_1012_count.csv', sep = ';')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()                                                                  #Step2
    parser.add_argument('-input1', type=str, required=True, help='gtf_file')  #Path of the input file   #Step3
    parser.add_argument('-output1', type=str, required=True, help='cassettes') #Path of the output file
    parser.add_argument('-output2', type=str, required=True, help='cassettes_count') 
    
    parser.add_argument('-input2', type=str, default = 'None', help='names of the nrRT sequences') 
    parser.add_argument('-output3', type=str, help='nr_cassettes') 
    parser.add_argument('-output4', type=str, help='nr_cassettes_count')

    args = parser.parse_args()                                                                          #Step4
    main(args)                                                                                          #Step5 (It is linked with 'def main(args)'. Could be: print('Hello,', args.name, args.lastname))
    
