import pandas as pd
from tqdm import tqdm
from absplice_scripts.data.DROP_annotations import sample_individual_mapping

df = pd.read_csv(snakemake.input['results'], sep='\t')

df_gene_map = pd.read_csv(snakemake.input['gene_map'], sep='\t')
gene_map = dict(zip(df_gene_map['gene_name'], df_gene_map['gene_id']))
df['gene_id'] = df['gene_name'].map(gene_map)

df = df[~df['gene_id'].isna()]

sample_map = pd.read_csv(snakemake.input['sample_map'], sep='\t')
sample_map = sample_individual_mapping(
    annotation_table_path=snakemake.input['sample_map'],
    key_assay=snakemake.params['key_assay'],                           
    value_assay=snakemake.params['value_assay'])

df['sample'] = df['Participant_ID'].map(sample_map)
df = df[~df['sample'].isna()]

df.to_csv(snakemake.output['results_annotated'], index=False)