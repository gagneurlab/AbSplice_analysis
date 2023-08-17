import pandas as pd
from absplice import SplicingOutlierResult

df_absplice_rna = pd.read_csv(snakemake.input['pred_absplice_rna'])
df_absplice_rna = df_absplice_rna.rename(columns={'y_pred': 'AbSplice_RNA'})

result = SplicingOutlierResult(
    df_absplice_rna = df_absplice_rna,
)

result.variant_absplice_rna.to_csv(snakemake.output['result_absplice_rna_variant'])
