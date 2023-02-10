import pandas as pd
from absplice import SplicingOutlierResult

result = SplicingOutlierResult(
    df_spliceai = snakemake.input['pred_spliceai'],
    gene_map = snakemake.input['gene_map'],
    df_var_samples = snakemake.input['var_samples_df']
)

result.gene_spliceai.to_parquet(snakemake.output['result_spliceai_gene'])

