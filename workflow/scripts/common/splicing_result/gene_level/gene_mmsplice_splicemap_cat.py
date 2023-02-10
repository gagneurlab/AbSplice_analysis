import pandas as pd
from absplice import SplicingOutlierResult

result = SplicingOutlierResult(
    df_mmsplice_cat = snakemake.input['pred_mmsplice_cat'],
)
result.gene_mmsplice_cat.to_parquet(snakemake.output['result_mmsplice_cat_gene'])

