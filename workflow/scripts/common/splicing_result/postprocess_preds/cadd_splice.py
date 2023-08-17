import pandas as pd
from absplice.utils import get_abs_max_rows

def variant_to_string(df):
    chrom = str(df['#Chrom'])
    if 'chr' not in chrom:
        chrom = 'chr' + chrom
    return chrom + ':' + str(df['Pos']) + ':' + df['Ref'] + '>' + df['Alt']

df_pred = pd.read_csv(snakemake.input['model'], sep='\t', skiprows=1)
df_pred = df_pred.rename(columns={'GeneID': 'gene_id', 'GeneName': 'gene_name'})
df_pred = df_pred[~df_pred['gene_id'].isna()]

df_pred['variant'] = df_pred.apply(lambda x: variant_to_string(x), axis=1)
df_pred = df_pred[['variant', 'gene_id', 'gene_name', 'RawScore', 'PHRED']]

# join sample info
if 'var_samples_df' in snakemake.input.keys():
    df_vcf_annotation = pd.read_csv(snakemake.input['var_samples_df'])
    
    df_pred = df_pred.set_index('variant').join(
        df_vcf_annotation.set_index('variant'), 
        how='inner')\
        .reset_index()

df_pred.to_parquet(
    snakemake.output['model_postprocess'], 
    index=False)