OUTPUT_DIR_RAW = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
OUTPUT_DIR_VAR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['variant_level']

import pandas as pd      
# --------------------------------variant level------------------------------- 
rule splicing_result_variant_level_absplice_dna_gtex_splicemaps:
    input:
        pred_absplice_dna = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['absplice']['absplice_preds']['gtex_splicemaps']['dna'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_absplice_dna_variant = OUTPUT_DIR_VAR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
    script:
        "../variant_level/variant_absplice_dna.py"
        
        
rule splicing_result_variant_level_absplice_rna_gtex_splicemaps:
    input:
        pred_absplice_rna = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['absplice']['absplice_preds']['gtex_splicemaps']['rna']['single_cat'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_absplice_rna_variant = OUTPUT_DIR_VAR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['single_cat'],
    script:
        "../variant_level/variant_absplice_rna.py"
        
        
rule splicing_result_variant_level_absplice_rna_gtex_splicemaps_all_cats:
    input:
        pred_absplice_rna = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['absplice']['absplice_preds']['gtex_splicemaps']['rna']['all_cats'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_absplice_rna_variant = OUTPUT_DIR_VAR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['all_cats'],
    script:
        "../variant_level/variant_absplice_rna.py"