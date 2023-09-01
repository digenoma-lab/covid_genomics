Preprocessing of COVID samples
================
2023-08-17

## Processing results in python

Cecret Nextflow pipeline was used to identify SARS-COV-2 variants in the
samples. This pipeline contains multiple tools that allows us to detect
variant types, mutations and heterogeneity in the samples. For this case
we use the outcome from freyja and iVar.

In this step we take the tsv files that contains the results of the
analysis of this tools and processing into tables to work with them
further.

``` python
import pandas as pd
import ast
import os 
#function to have the full path of each file in a folder
def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]
  
#position variant from freya
freya_files = sorted(listdir_fullpath('join_flow_cells/freyja'))

with open('freya_results.txt', 'w') as f: #file to save results
    f.write('sampleId,th_filter,pass_filter')
    f.write('\n')
    for file in freya_files:
        
        if file.endswith('variants.tsv'): #open files where are the results of variants per position 
            data = pd.read_csv(file, sep='\t')
            #choosing to type of filters to select data
            th_filter = len(data[data["ALT_DP"]>400])
            pass_filter = len(data[data["PASS"] == True])
            sample_name = file.split('.')[0].split('/')[-1:][0][:-9] #sample name 
            f.write('{},{},{}'.format(sample_name,th_filter,pass_filter))
            f.write('\n')
```

``` python

#covid variant 
covid_variant = pd.read_csv('join_flow_cells/freyja/aggregated-freyja.tsv', sep='\t')
all_variants = {}
for index, row in covid_variant.iterrows():
    sample = 'sample_{}'.format(row['Unnamed: 0'].split('.')[0].split('_')[2])
    variants = ast.literal_eval(row['summarized'])
    all_variants[sample] = dict(variants)

all_variants = pd.DataFrame(all_variants)
all_variants.to_csv('freya_variants.tsv', sep='\t')

ivar_files = sorted(listdir_fullpath('join_flow_cells/ivar_variants'))
positions = {}

for file in ivar_files:   
    if file.endswith('variants.tsv'): #open files where are the results of variants per position 
        sample_name = 'sample_{}'.format(file.split('.')[0].split('/')[-1][-2:])
        data = pd.read_csv(file, sep='\t')
        count = {}
        for pos in data[['POS']].values:
            count[pos[0]] = 1
       
        positions[sample_name] = count

ivar_results = pd.DataFrame(positions).fillna(0).transpose()

ivar_results.to_csv('ivar_results.tsv', sep='\t', float_format='%.0f')

#
samtools_depth = sorted(listdir_fullpath('join_flow_cells/samtools_depth'))

samtools_data = {}
for file in samtools_depth:      
    if file.endswith('depth.txt'):
        file_name = file.split('.')[0].split('/')[-1:][0].split('_')[-1]
        sample_name = 'sample_{}'.format(file_name)[:-3]
        data = pd.read_csv(file, sep='\t', names=['ref', 'pos', 'depth'])
        #transforming data into a dictionary for easy transform into dataframe
        values = pd.Series(data.depth.values,index=data.pos).to_dict()
        samtools_data[sample_name] = values
        
pd.DataFrame(samtools_data).to_csv('samtools_depth.csv', sep=',')
```

## Joining data in R

After python codes generate tables from the iVar and Freyja output this
are joined using R with the demographic information.

    ## Using X as id variables

    ##    sampleId pass_freya mut_ivar variant location age analisis_data
    ## 1 sample_01       6330       15   Other Rancagua  28    2020-07-17
    ## 2 sample_02       8613       15   Other Rancagua  32    2020-08-04
    ## 3 sample_03       2618       24   Other Rancagua  44    2020-08-14
    ## 4 sample_04      25573       15   Other Rancagua  24    2020-08-19
    ## 5 sample_05      22013       12   Other Rancagua  56    2020-08-25
    ## 6 sample_06      16678       17   Other Rancagua  36    2020-09-02
