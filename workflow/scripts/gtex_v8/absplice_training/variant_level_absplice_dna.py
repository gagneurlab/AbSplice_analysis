import pandas as pd
from absplice import SplicingOutlierResult

df_absplice_dna = pd.read_csv(snakemake.input['pred_absplice_dna'])
df_absplice_dna = df_absplice_dna.rename(columns={'y_pred': 'AbSplice_DNA'})

result = SplicingOutlierResult(
    df_absplice_dna = df_absplice_dna,
)

result.variant_absplice_dna.to_csv(snakemake.output['result_absplice_dna_variant'])
