OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
ABSPLICE_INPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['combined_input']

import pandas as pd        
# --------------------------------predict AbSplice------------------------------- 
rule splicing_pred_absplice_dna_gtex_splicemaps:
    input:
        absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
        absplice_model = config_precomputed['absplice']['dna'],
    params:
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        absplice_dna_pred = OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
    script:
        "../absplice_dna.py"
        
        
rule splicing_pred_absplice_rna_gtex_splicemaps:
    input:
        absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['single_cat'],
        absplice_model = config_precomputed['absplice']['rna'],
    params:
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        absplice_rna_pred = OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['single_cat'],
    script:
        "../absplice_rna.py"
        
        
rule splicing_pred_absplice_rna_gtex_splicemaps_all_cats:
    input:
        absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['all_cats'],
        absplice_model = config_precomputed['absplice']['rna'],
    params:
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        absplice_rna_pred = OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['all_cats'],
    script:
        "../absplice_rna.py"