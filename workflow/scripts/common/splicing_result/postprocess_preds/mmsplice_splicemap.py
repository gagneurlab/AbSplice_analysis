import pandas as pd
from absplice import SplicingOutlierResult

if 'var_samples_df' in snakemake.input.keys():
    result = SplicingOutlierResult(
        df_mmsplice = snakemake.input['model'],
        gene_tpm = snakemake.input['gene_tpm'],
        df_var_samples = snakemake.input['var_samples_df'],
    )
else:
    result = SplicingOutlierResult(
        df_mmsplice = snakemake.input['model'],
        gene_tpm = snakemake.input['gene_tpm'],
    )

result.df_mmsplice.to_parquet(
    snakemake.output['model_postprocess'], 
    index=False)

