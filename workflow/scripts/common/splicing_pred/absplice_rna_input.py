import pandas as pd
from absplice import SplicingOutlierResult
from absplice.utils import get_abs_max_rows

groupby_index = ['variant', 'gene_id', 'sample', 'tissue_cat']

df_outliers_cat = pd.concat([pd.read_csv(i) for i in snakemake.input['CAT_pval']])
df_outliers_cat = get_abs_max_rows(
    df_outliers_cat.set_index(groupby_index), 
    groupby_index, 
    'pValueGene_g_minus_log10'
).reset_index()

result = SplicingOutlierResult(
    df_absplice_dna_input = snakemake.input['absplice_dna_input'],
    df_mmsplice_cat = snakemake.input['pred_mmsplice_cat'],
    df_outliers_cat = df_outliers_cat,
)

result.absplice_rna_input.to_parquet(snakemake.output['absplice_rna_input'])
