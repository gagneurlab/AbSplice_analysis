import pandas as pd
from absplice import SplicingOutlierResult

result = SplicingOutlierResult(
    df_absplice_rna = snakemake.input['pred_absplice_rna'],
)

result.gene_absplice_rna.to_parquet(snakemake.output['result_absplice_rna_gene'])
