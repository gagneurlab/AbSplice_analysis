import pandas as pd
from absplice import SplicingOutlierResult

result = SplicingOutlierResult(
    df_mmsplice = snakemake.input['pred_mmsplice'],
    df_spliceai = snakemake.input['pred_spliceai'],
    gene_map = snakemake.input['gene_map'],
    gene_tpm = snakemake.input['gene_tpm'],
    df_var_samples = snakemake.input['var_samples_df']
)

result.absplice_dna_input.to_parquet(snakemake.output['absplice_dna_input'])
