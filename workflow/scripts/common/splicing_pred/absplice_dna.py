import pandas as pd
from absplice import SplicingOutlierResult

features_dna = sorted([
    'delta_psi', 'delta_score', 'splice_site_is_expressed', 'delta_logit_psi'
])

result = SplicingOutlierResult(
    df_absplice_dna_input = snakemake.input['absplice_dna_input'],
)

result.predict_absplice_dna(
    pickle_file=snakemake.input['absplice_model'],
    median_n_cutoff=snakemake.params['median_n_cutoff'],  
    # tpm_cutoff=snakemake.params['tpm_cutoff'],
    features=features_dna, 
    abs_features=False,
    extra_info=False    
)

# result._absplice_dna.to_csv(snakemake.output['absplice_dna_pred'])
result._absplice_dna.to_parquet(snakemake.output['absplice_dna_pred'], engine='pyarrow')
