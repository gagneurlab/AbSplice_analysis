import pandas as pd
from absplice import SplicingOutlierResult
from absplice.utils import get_abs_max_rows
from absplice_scripts.utils.mapping_utils import subset_tissues

result = SplicingOutlierResult(
    df_mmsplice = snakemake.input['pred_mmsplice'],
    df_var_samples = snakemake.input['var_samples_df'],
    gene_tpm = snakemake.input['gene_tpm']
)

df = result.df_mmsplice

tissue_map = pd.read_csv(snakemake.input['tissue_map'])
tissue_map = dict(zip(tissue_map['tissue'], tissue_map['tissue_main']))
df = subset_tissues(df, tissue_map, chosen_tissue=snakemake.wildcards['gtex_tissue'])

df = get_abs_max_rows(
    df = df.set_index(['gene_id', 'sample']),
    groupby = ['gene_id', 'sample'],
    max_col = 'delta_logit_psi',
)

assert df.shape[0] == len(df.index.unique())

df.to_parquet(snakemake.output['result_mmsplice_gene'])

