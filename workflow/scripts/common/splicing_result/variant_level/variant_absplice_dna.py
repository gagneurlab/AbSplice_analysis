import pandas as pd
from absplice import SplicingOutlierResult
from absplice.utils import get_abs_max_rows
from absplice_scripts.utils.mapping_utils import subset_tissues

df = pd.read_csv(snakemake.input['pred_absplice_dna'])

df = get_abs_max_rows(
    df = df.set_index(['variant', 'gene_id', 'sample', 'tissue']),
    groupby = ['variant', 'gene_id', 'sample', 'tissue'],
    max_col = 'AbSplice_DNA',
#     dropna=False
)

assert df.shape[0] == len(df.index.unique())

df.to_csv(snakemake.output['result_absplice_dna_variant'])

