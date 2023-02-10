import pandas as pd
from tqdm import tqdm
import pandas as pd
from absplice.utils import get_abs_max_rows

def variant_to_string(df):
    chrom = str(df['#Chrom'])
    if 'chr' not in chrom:
        chrom = 'chr' + chrom
    return chrom + ':' + str(df['Pos']) + ':' + df['Ref'] + '>' + df['Alt']

df_pred = pd.read_csv(snakemake.input['pred_cadd_splice'], sep='\t', skiprows=1)
df_pred = df_pred.rename(columns={'GeneID': 'gene_id', 'GeneName': 'gene_name'})
df_pred = df_pred[~df_pred['gene_id'].isna()]
df_vcf_annotation = pd.read_csv(snakemake.input['var_samples_df'])

df_pred['variant'] = df_pred.apply(lambda x: variant_to_string(x), axis=1)
df_pred = df_pred[['variant', 'gene_id', 'gene_name', 'RawScore', 'PHRED']]

# join sample info
df_pred = df_pred.set_index('variant').join(
    df_vcf_annotation.set_index('variant'), 
    how='inner')\
    .reset_index()

# Aggregate predictions
df_pred_agg = get_abs_max_rows(
    df_pred.set_index(['variant', 'gene_id', 'sample']), 
    groupby=['variant', 'gene_id', 'sample'], 
    max_col=snakemake.params['max_col']
)

df_pred_agg.to_csv(snakemake.output['result_cadd_splice_variant'])
