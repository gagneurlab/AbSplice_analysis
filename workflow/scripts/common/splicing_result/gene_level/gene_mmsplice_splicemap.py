import pandas as pd
from absplice import SplicingOutlierResult
from absplice.utils import get_abs_max_rows

result = SplicingOutlierResult(
    df_mmsplice = snakemake.input['pred_mmsplice'],
    df_var_samples = snakemake.input['var_samples_df'],
    gene_tpm = snakemake.input['gene_tpm']
)

# Aggregate predictions
df_pred_agg = get_abs_max_rows(
    result.df_mmsplice.set_index(['gene_id', 'sample', 'tissue']), 
    groupby=['gene_id', 'sample', 'tissue'], 
    max_col='delta_logit_psi'
)

df_pred_agg.to_parquet(snakemake.output['result_mmsplice_gene'])

