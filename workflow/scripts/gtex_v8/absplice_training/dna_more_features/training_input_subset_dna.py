import pandas as pd
import numpy as np

df = pd.read_csv(snakemake.input['absplice_input'])
df = df[
    (~df['gene_id'].isna())
    & (
        (np.abs(df['delta_psi']) > 0.01)
        | (df['delta_score'] > 0.01)
        | (df['PHRED'] > 5)
        | (df['squirls_scores'] > 0.01)
        | (df['delta_logit_psi_mtsplice'] > 1)
    )
]
df.to_csv(snakemake.output['absplice_input'], index=False)
