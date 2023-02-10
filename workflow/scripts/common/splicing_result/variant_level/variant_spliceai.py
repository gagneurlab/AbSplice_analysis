import pandas as pd
from absplice import SplicingOutlierResult
from absplice.utils import get_abs_max_rows

result = SplicingOutlierResult(
    df_spliceai = snakemake.input['pred_spliceai'],
    gene_map = snakemake.input['gene_map'],
    df_var_samples = snakemake.input['var_samples_df']
)

df = result.df_spliceai

df = get_abs_max_rows(
    df = df.set_index(['variant', 'gene_id', 'sample']),
    groupby = ['variant', 'gene_id', 'sample'],
    max_col = 'delta_score',
#     dropna=False
)

assert df.shape[0] == len(df.index.unique())

df.to_csv(snakemake.output['result_spliceai_variant'])

