import pandas as pd
from absplice import SplicingOutlierResult
from absplice.utils import get_abs_max_rows
from absplice_scripts.utils.mapping_utils import subset_tissues

result = SplicingOutlierResult(
    df_mmsplice_cat = snakemake.input['pred_mmsplice_cat'],
)

df = result.df_mmsplice_cat

df = df[~df['delta_psi_cat'].isna()]

tissue_map = pd.read_csv(snakemake.input['tissue_map'])
tissue_map = dict(zip(tissue_map['tissue'], tissue_map['tissue_main']))
df = subset_tissues(df, tissue_map, chosen_tissue=snakemake.params['gtex_tissue'])

df = get_abs_max_rows(
    df = df.set_index(['variant', 'gene_id', 'sample']),
    groupby = ['variant', 'gene_id', 'sample'],
    max_col = 'delta_psi_cat',
)

assert df.shape[0] == len(df.index.unique())

df.to_csv(snakemake.output['result_mmsplice_cat_variant'])

