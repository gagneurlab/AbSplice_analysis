import pandas as pd
import numpy as np
from absplice.utils import get_abs_max_rows

def minus_log10(x):
    return -np.log10(x)

df_outliers = pd.read_csv(snakemake.input['outliers_variant_level'])
df_outliers['tissue_cat'] = snakemake.wildcards['tissue']

df = pd.read_csv(snakemake.input['outliers_pval'])
df['pValueGene_g_minus_log10'] = df['pValueGene_g'].apply(lambda x: minus_log10(x))
df = get_abs_max_rows(
    df.set_index(['variant', 'gene_id', 'sample']), 
    ['variant', 'gene_id', 'sample'], 
    'pValueGene_g_minus_log10'
).reset_index()

join_index = ['variant', 'gene_id', 'sample']
df_outliers = df_outliers.set_index(join_index).join(
    df.set_index(join_index)).reset_index()

df_outliers = df_outliers[[
    'variant', 'gene_id', 'sample',
    'pValueGene_g_minus_log10', 'tissue_cat'
]]

df_outliers.to_csv(snakemake.output['outliers_pval'], index=False)
    