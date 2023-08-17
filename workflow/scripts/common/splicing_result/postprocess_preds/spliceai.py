import pandas as pd
from absplice import SplicingOutlierResult

if 'var_samples_df' in snakemake.input.keys():
    result = SplicingOutlierResult(
        df_spliceai = snakemake.input['model'],
        gene_map = snakemake.input['gene_map'],
        df_var_samples = snakemake.input['var_samples_df']
    )
else:
    result = SplicingOutlierResult(
        df_spliceai = snakemake.input['model'],
        gene_map = snakemake.input['gene_map'],
    )

result.df_spliceai.to_parquet(
    snakemake.output['model_postprocess'],
    index=False)

