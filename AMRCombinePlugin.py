# The purpose of the script is the following:
#   1) combine the AMR gene counts into a single file;
#   2) Count the number of unique genes in each cohort;
#   3) Plot the counts of AMR genes across multiple time points;

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

import PyPluMA
import PyIO

##################################### INPUT #######################################
class AMRCombinePlugin:
  def input(self, inputfile):
    self.parameters = PyIO.readParameters(inputfile)
  def run(self):
    pass
  def output(self, outputfile):
   samples_list = [x.strip('\n') for x in open(PyPluMA.prefix()+"/"+self.parameters["samples"]).readlines()]
   results_dir = PyPluMA.prefix()+"/"+self.parameters["results"]+"/"
   metadata_file = PyPluMA.prefix()+"/"+self.parameters["metadata"]
###################################################################################

   print("Combining ABR abundance for {} samples...".format(len(samples_list)))

   # The final table with all ABR genes combined:
   combined_df = pd.DataFrame({'sampleID':[], 'gene_name':[], 'ARG':[], 'read_count':[], 'gene_length':[],'coverage':[]})

   for sample_id in samples_list:
    result_i = results_dir+sample_id + '/' + 'groot_out_' + sample_id +'.txt'
    df_i = pd.read_csv(result_i, sep='\t', names=['ARG', 'read_count', 'gene_length','coverage'])
    df_i['gene_name'] = df_i['ARG'].apply(lambda x: x.split('|')[-1] if x.split('|')[-1]!='RequiresSNPConfirmation' else x.split('|')[-2])

    # remove duplicates by selecting genes with maximum read counts:
    unique_df_i = df_i[df_i.groupby(['gene_name'])['read_count'].transform(max) == df_i['read_count']]

    # sometimes the read_count is exactly the same for several gene variants. In this case, randomly remove the duplicates:
    unique_df_i = unique_df_i.drop_duplicates(subset=['gene_name'])

    # add column with sample IDs
    unique_df_i['sampleID'] = sample_id
    unique_df_i = unique_df_i[['sampleID', 'gene_name', 'ARG', 'read_count', 'gene_length','coverage']]
    combined_df = combined_df.append(unique_df_i)

   combined_df = combined_df.reset_index(drop=True)
   #combined_df["sampleID"] = combined_df["sampleID"].apply(lambda x: x[:4])

   print("Total unique ABR genes identified: {}".format(len(combined_df['gene_name'].unique())))

   # Combining the file with metdata:
   metadata_df = pd.read_csv(metadata_file)

   combined_metadata_df = combined_df.merge(metadata_df, how='left', left_on='sampleID', right_on="Samples#")

   #if not os.path.exists('amr_counts_metagen'):
   #os.mkdir('amr_counts_metagen')

   combined_metadata_df.to_csv(outputfile, index=False)
