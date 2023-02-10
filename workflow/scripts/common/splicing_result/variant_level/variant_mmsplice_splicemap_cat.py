import pandas as pd
from absplice import SplicingOutlierResult
from absplice.utils import get_abs_max_rows
from absplice_scripts.utils.mapping_utils import subset_tissues

result = SplicingOutlierResult(
    df_mmsplice_cat = snakemake.input['pred_mmsplice_cat'],
)

df = result.df_mmsplice_cat

df = df[~df['delta_psi_cat'].isna()]

df = get_abs_max_rows(
    df = df.set_index(['variant', 'gene_id', 'sample', 'tissue']),
    groupby = ['variant', 'gene_id', 'sample', 'tissue'],
    max_col = 'delta_psi_cat',
#     dropna = False,
)

assert df.shape[0] == len(df.index.unique())

df.to_csv(snakemake.output['result_mmsplice_cat_variant'])

