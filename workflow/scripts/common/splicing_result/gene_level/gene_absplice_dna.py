import pandas as pd
from absplice import SplicingOutlierResult

result = SplicingOutlierResult(
    df_absplice_dna = snakemake.input['pred_absplice_dna'],
)

result.gene_absplice_dna.to_parquet(snakemake.output['result_absplice_dna_gene'])
