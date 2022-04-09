"""
Python Script
For: Create table with VRs and ORFs, add column 'VR_position', and then add the Biome and Taxonomy information to each target gene.

Inputs & Outputs:
- input1 --> './test6_magswithDGRgtf/'    
- input2 --> './test6_prokkaoutputs'
- output1 --> './test6_funannotation.csv'

- input3 --> './test6_outputfiles/nayfach2020_genome_metadata_curated_biomtaxa.tsv'
- output2 --> './test6_funannotation_complete_taxa.csv'

# Example of how to run the script:
#python .\script_functional_annotation.py -input1 '' -input2 '' -output1 '' -input3 '' -output2 ''

"""

import argparse             #Step1
import os 
from os import listdir
from os.path import isfile, join
import pandas as pd
from tqdm import tqdm


def main(args):
    path_gtfVR = args.input1
    path_gffORF = args.input2
    path_output1 = args.output1
    path_metadata = args.input3
    path_output2 = args.output2
    

    #####   1. From each _out-DGR.gtf in ./test6_magswithDGRgtf, extract the information of the VR columns.
    #####   Extract the information of the columns (in .gtf): #genome_id, #contig_id, #VR_strand, #VR_start, #VR_end

    # Set my directory with input files (.gtf).
    path_VR = path_gtfVR                                                    #1 ARGS   #path_gtfVR = './test6_magswithDGRgtf/'
    files_VR = [join(path_VR, f) for f in listdir(path_VR) if isfile(join(path_VR, f))]

    # Create an empty list for storing the .gtf dataframes.
    df_VR = []

    # For each .gtf extract the #genome_id, and add it to a new column on each .gtf.
    for path_file in files_VR:
        genomeid = path_file.split('/')     #['.','test6_magswithDGR','33000000_001_out-DGR.gtf']
        genomeid = genomeid[-1]             #'33000000_001_out-DGR.gtf'       
        genomeid = genomeid[:-12]           #'33000000_001'
        df = pd.read_csv(path_file, sep = '\t', header=None, encoding='latin', engine = 'python')
        df['genome_id'] = genomeid

    # Here I will append my .gtf dataframes to the "df_VR" list.
        df_VR.append(df)

    # Here I will join all my .gtf dataframes from my "df_VR" list in a single dataframe "df_VR".
    df_VR = pd.concat(df_VR)
    df_VR = df_VR.reset_index()
    #print(df_VR)

    asterisk_items = df_VR[df_VR[3] == '*']

    asterisk_items = asterisk_items.index.tolist()

    for index in asterisk_items:
        value_start = df_VR.iloc[index - 1][3]
        value_end = df_VR.iloc[index - 1][4]
        df_VR.loc[index, 3] = value_start
        df_VR.loc[index, 4] = value_end
    #print(df_VR)

    # Here I will filtrate the columns I want: #genome_id, #contig_id, #VR, #VR_strand, #TR_orig_start, #TR_orig_end, #VR_start, #VR_end, #A-N mut, #non A-N mut
    # I will keep the rows with the value == VR on the 1st column.
    # I will rename the columns and change the numbers to integers.
    cols = ['genome_id',0,1,2,3,4,5,6,7,8]
    df_VR = df_VR[cols]
    df_VR = df_VR[df_VR[1] == 'VR']
    df_VR.columns = ['genome_id', 'contig_id', 'VR', 'VR_strand', 'TR_orig_start', 'TR_orig_end', 'VR_start', 'VR_end', 'A-N mut', 'non A-N mut']
    df_VR = df_VR.astype({'TR_orig_start': int, 'TR_orig_end': int, 'VR_start': int, 'VR_end': int, 'A-N mut' : int, 'non A-N mut': int})
    #print(df_VR)

    # Here I will make a list with the values of the #genome_id column from my "df_VR" dataframe.
    genomeids = list(df_VR['genome_id'].unique())

    #####   2. In the directory ./test6_prokkaoutputs, search the #genome_id.gff and make column #ORF_id, also search the #genomeid.tsv
    #####   Extract the information of the columns (in .gff): #contig_id, #ORF_strand, #ORF_start, #ORF_end, #ORF_id

    # Set my directory with input files (.gff y .tsv)
    path_ORF = path_gffORF                                                                #2 ARGS   #path_gffORF = './test6_prokkaoutputs'

    df_gff = []
    df_tsv =[]
    for genomeid in tqdm(genomeids):
        path_gff = f'{path_ORF}/{genomeid}/{genomeid}.gff'
        if isfile(path_gff):
            df = pd.read_csv(path_gff, sep = '\t', header = None, comment = '#', encoding='latin', engine = 'python')  #comment = "#" va a omitir las filas que empiecen con #
            df_gff.append(df)

        path_tsv = f'{path_ORF}/{genomeid}/{genomeid}.tsv'
        if isfile(path_tsv):
            df = pd.read_csv(path_tsv, sep = '\t', encoding='latin', engine = 'python')
            df_tsv.append(df)
        
    df_gff = pd.concat(df_gff)
    df_tsv = pd.concat(df_tsv)
    #print(df_gff)
    #print(df_tsv)

    # Create a function for obtaining the ORF_id from one of the colummns in the .gff dataframe.
    def get_id(rawcol):
        rawcol = rawcol.split(';')  #ID=3300027819_166_00001;prediction;product
        id = rawcol[0] #ID=3300027819_166_00001
        id = id[3:] #3300027819_166_00001
        return id

    # Here I will filtrate the columns I want: 'contig_id', 'database', 'ORF_ftype', 'ORF_start', 'ORF_end', 'ORF_strand', 'ORF_id'.
    # I will apply the function 'get_id' for getting the 'ORF_id' on that column.
    cols = [0,1,2,3,4,6,8]
    df_gff = df_gff[cols]
    df_gff = df_gff.dropna()
    df_gff[8] = df_gff[8].apply(get_id)
    df_gff.columns = ['contig_id', 'database', 'ORF_ftype', 'ORF_start', 'ORF_end', 'ORF_strand', 'ORF_id']
    #print(df_gff)

    #print('archivo gff:', df_gff.columns)
    #print('archivo VR:', df_VR.columns)

    #####   3. Search the column #contig_id in my files .gtf (VRs) and .gff (ORFs) and merge them.

    # Add the information from the .gtf (VR) dataframe, using the column 'contig_id' as the link.
    df_gff = pd.merge(df_gff, df_VR, on = 'contig_id', how = 'left') #how='left' will keep all the rows of the "df_gff", and duplicate some of them if necessary.
    df_gff = df_gff.dropna() #this will drop 'contig_id' rows where there is no information of VRs.
    #print(df_gff[['contig_id', 'VR_start', 'VR_end', 'ORF_start', 'ORF_end', 'ORF_id']])

    #####   4. Check if the #VR_start or #VR_end is between the #ORF_start and #ORF_end.

    # Once I have the 'contig_id' rows that have VR information, I need to check if the VR is within the ORF.
    df_gff = df_gff[((df_gff['ORF_start'] < df_gff['VR_start']) & (df_gff['VR_start'] < df_gff['ORF_end'])) |
                    ((df_gff['ORF_start'] < df_gff['VR_end']) & (df_gff['VR_end'] < df_gff['ORF_end']))]

    #df_gff = df_gff[(df_gff['ORF_start'] < df_gff['VR_start']) & (df_gff['VR_end'] < df_gff['ORF_end'])]

    #print(df_gff[['contig_id', 'VR_start', 'VR_end', 'ORF_start', 'ORF_end', 'ORF_id']])
    #print(df_tsv.columns)
    #print(df_gff.columns)

    #####   5. In the directory ./test6_prokkaoutputs, search the #genomeid.tsv.
    #####   Search the column #ORF_id in my files .gff and .tsv and merge them.
    #####   Extract the information of the columns (in .tsv): #ORF_id #ORF_length, #ORF_ftype, #ORF_gene, #ORF_COG, #ORF_product

    # Add the information from the .tsv dataframe, using the column 'ORF_id' as the link.
    df_tsv = df_tsv.rename({'locus_tag': 'ORF_id'}, axis = 1)
    df_gff = pd.merge(df_gff, df_tsv, on = 'ORF_id', how = 'left')

    # Save the new table with the VRs (end, start) and the ORFs (end, start, other info)
    df_gff.rename({'length_bp': 'ORF_length_bp', 'A-N mut': 'VR_A-N mut', 'non A-N mut':'VR_non A-N mut'}, axis = 'columns', inplace=True)
    cols = ['genome_id', 'contig_id', 'ORF_id', 'ORF_length_bp', 'ORF_strand', 'ORF_start', 'ORF_end', 'VR', 'VR_strand', 'VR_start', 'VR_end', 'TR_orig_start', 'TR_orig_end', 'VR_A-N mut', 'VR_non A-N mut', 'database', 'ORF_ftype', 'ftype', 'gene', 'EC_number', 'COG', 'product']
    df_gff = df_gff[cols]
    #print(df_gff.columns)

    #####   6. Check if columns 'ORF_strand' and 'VR_strand' are equal ( + == +, - == -). Keep only these rows.
    #####   https://stackoverflow.com/questions/43951558/remove-rows-that-two-columns-have-the-same-values-by-pandas

    df_gff = df_gff[df_gff['ORF_strand'] == df_gff['VR_strand']]

    #####   7. Add columns 'VR_length_bp' and 'TR_length_bp'
    VR_length_bp = df_gff['VR_end'] - df_gff['VR_start'] + 1
    TR_length_bp = df_gff['TR_orig_end'] - df_gff['TR_orig_start'] + 1

    df_gff['VR_length_bp'] = VR_length_bp
    df_gff['TR_length_bp'] = TR_length_bp

    cols2 = ['genome_id', 'contig_id', 'ORF_id', 'ORF_length_bp', 'ORF_strand', 'ORF_start', 'ORF_end', 'VR', 'VR_strand', 'VR_start', 'VR_end', 'TR_orig_start', 'TR_orig_end', 'VR_length_bp', 'TR_length_bp', 'VR_A-N mut', 'VR_non A-N mut', 'database', 'ORF_ftype', 'ftype', 'gene', 'EC_number', 'COG', 'product']
    df_gff = df_gff[cols2]

    #####   8. Add column 'VR_position'
    #ORF_start, ORF_end --> 3 sections (begin, middle, end)
    #If ORF_strand == + --> position of VR_start
    #If ORF_strand == - --> position of VR_end
    #Add new column 'VR_position' (3 possibilities)

    #print(df_gff[['ORF_length_bp', 'ORF_strand', 'ORF_start', 'ORF_end', 'VR_start', 'VR_end']])

    VR_position = []

    for index, row in df_gff.iterrows():
        ORF_start = row['ORF_start']
        ORF_end = row['ORF_end']
        ORF_length = row['ORF_length_bp']
        part_length = ORF_length/3
        cut1 = ORF_start + part_length
        cut2 = cut1 + part_length
        ORF_strand = row['ORF_strand']
        value = None
        if ORF_strand == '-':
            value = row['VR_end']
        else:
            value = row['VR_start']

        position = None
        if ORF_start < value <= cut1:
            position = 'begin'
        if cut1 < value <= cut2:
            position = 'middle'
        if cut2 < value < ORF_end:
            position = 'end'
        
        VR_position.append(position)

    df_gff['VR_position'] = VR_position

    #print(df_gff[['ORF_length_bp', 'ORF_strand', 'ORF_start', 'ORF_end', 'VR_start', 'VR_end', 'VR_position']])

    #df_negative = df_gff[df_gff['ORF_strand'] == '-']
    #df_positive = df_gff[df_gff['ORF_strand'] == '+']

    cols3 = ['genome_id', 'contig_id', 'ORF_id', 'ORF_length_bp', 'ORF_strand', 'ORF_start', 'ORF_end', 'VR', 'VR_strand', 'VR_start', 'VR_end', 'VR_position', 'TR_orig_start', 'TR_orig_end', 'VR_length_bp', 'TR_length_bp', 'VR_A-N mut', 'VR_non A-N mut', 'database', 'ORF_ftype', 'ftype', 'gene', 'EC_number', 'COG', 'product']
    df_gff = df_gff[cols3]

    #df_gff.to_csv('./test6_funannotation.csv', index = False)                            #3 ARGS
    df_gff.to_csv(path_output1, index = False)                                            #3 ARGS   path_output1 = './test6_funannotation.csv'

    #============================================================================================
    # I want to add information from .tsv "nayfach2020_genome_metadata_curated_biomtaxa.tsv"
    # I want to add information from columns 'biome', 'ecosystem_category', 'ecosystem_type', 'habitat'
    
    if args.input3 != 'None':    #if the file exists (not None), run this part of the script
        path_fulldata = path_metadata                                                         #4 ARGS   #path_metadata = './test6_outputfiles/nayfach2020_genome_metadata_curated_biomtaxa.tsv'
        df_fulldata = pd.read_csv(path_fulldata, sep='\t')

        # Keep only some columns
        cols = ['genome_id','metagenome_id','genome_length','num_contigs','otu_id','biome','ecosystem_category','ecosystem_type', 'habitat', 'domain',  'phylum']
        df_fulldata = df_fulldata[cols]
        #print(df_fulldata) 

        # Join dataframes, using the column 'genome_id'
        df_gff = pd.merge(df_gff, df_fulldata, on = 'genome_id', how = 'left') #how='left' will keep all the rows of the "df_gff", and duplicate some of them if necessary.
        #df_gff = df_gff.dropna() 

        # Save the new table VRs/ORFs/BIOMES 
        #df_gff.to_csv('./test6_funannotation_complete_taxa.csv', index = False)                            #5 ARGS     
        df_gff.to_csv(path_output2, index = False)                                                          #5 ARGS  #path_output2 = './test6_funannotation_complete_taxa.csv'



if __name__ == '__main__':
    parser = argparse.ArgumentParser()                                                                               #Step2
    parser.add_argument('-input1', type=str, required=True, help='directory_gtfDGRs')  #Path of the input file       #Step3
    parser.add_argument('-input2', type=str, required=True, help='directory_prokkaoutputs') 
    parser.add_argument('-output1', type=str, required=True, help='functional_annotation.csv') #Path of the output file
    
    parser.add_argument('-input3', type=str, default = 'None', help='metadata.tsv') 
    parser.add_argument('-output2', type=str, help='functional_annotation_complete.csv') 

    args = parser.parse_args()                                                                          #Step4
    main(args)                                                                                          #Step5 (It is linked with 'def main(args)')
