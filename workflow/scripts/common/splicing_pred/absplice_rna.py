import pandas as pd
from absplice import SplicingOutlierResult

features_rna = sorted([
    'delta_psi', 'delta_score', 'splice_site_is_expressed', 'delta_logit_psi', 
    'delta_psi_cat', 'pValueGene_g_minus_log10'
])
    
result = SplicingOutlierResult(
    df_absplice_rna_input = snakemake.input['absplice_rna_input'],
)

result.predict_absplice_rna(
    pickle_file=snakemake.input['absplice_model'],
    median_n_cutoff=snakemake.params['median_n_cutoff'],  
    tpm_cutoff=snakemake.params['tpm_cutoff'],
    features=features_rna, 
    abs_features=False
)

result._absplice_rna.to_parquet(snakemake.output['absplice_rna_pred'])
