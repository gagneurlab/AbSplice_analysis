import pandas as pd
from absplice import SplicingOutlierResult
from absplice.utils import get_abs_max_rows
from absplice_scripts.utils.mapping_utils import subset_tissues

result = SplicingOutlierResult(
    df_mmsplice = snakemake.input['pred_mmsplice'],
    df_var_samples = snakemake.input['var_samples_df']
)

df = result.df_mmsplice

df = get_abs_max_rows(
    df = df.set_index(['variant', 'gene_id', 'sample', 'tissue']),
    groupby = ['variant', 'gene_id', 'sample', 'tissue'],
    max_col = 'delta_psi',
)

assert df.shape[0] == len(df.index.unique())

df.to_parquet(snakemake.output['result_mmsplice_variant'])

