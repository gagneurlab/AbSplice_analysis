OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
ABSPLICE_INPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['combined_input']

import pandas as pd
# --------------------------------AbSplice input-------------------------------         
rule splicing_pred_absplice_dna_input_gtex_splicemaps_no_samples:
    input:
        pred_mmsplice = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap']['gtex_splicemaps'],
        pred_spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
        gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
        gene_tpm = config_precomputed['gene_tpm_gtex'].format(gtex_version=config['gtex_version']),
        splicemap_5 = expand(config_precomputed['splicemap_gtex']['psi5'],
                             gtex_version=config['gtex_version'], tissue=config['gtex_tissues']),
        splicemap_3 = expand(config_precomputed['splicemap_gtex']['psi3'],
                             gtex_version=config['gtex_version'], tissue=config['gtex_tissues']),
    params:
        gtex_tissues = config['gtex_tissues']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 48000,
        threads = 1
    output:
        absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
    script:
        "../absplice_dna_input_no_samples.py"
        
        
# --------------------------------predict AbSplice------------------------------- 
rule splicing_pred_absplice_dna_gtex_splicemaps:
    input:
        absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
        absplice_model = config_precomputed['absplice']['dna'],
    params:
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 48000,
        threads = 1
    output:
        absplice_dna_pred = OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
    script:
        "../predict_absplice_dna.py"