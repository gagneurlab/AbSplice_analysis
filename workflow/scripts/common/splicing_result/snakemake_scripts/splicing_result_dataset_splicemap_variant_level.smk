# --------------------------------variant level------------------------------- 
rule splicing_result_variant_level_mmsplice_dataset_splicemap:
    input:
        pred_mmsplice = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['raw_pred']['mmsplice_splicemap']['dataset_splicemap'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        result_mmsplice_variant = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['pred_variant_level']['mmsplice_splicemap']['dataset_splicemap']
    script:
        "../variant_level/variant_mmsplice_splicemap.py"
        
        
rule splicing_result_variant_level_mmsplice_cat_dataset_splicemap:
    input:
        pred_mmsplice_cat = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['raw_pred']['mmsplice_splicemap_cat']['dataset_splicemap'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        result_mmsplice_cat_variant = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['pred_variant_level']['mmsplice_splicemap_cat']['dataset_splicemap'],
    script:
        "../variant_level/variant_mmsplice_splicemap_cat.py"
        
       
    
rule splicing_result_variant_level_absplice_dna_dataset_splicemap:
    input:
        pred_absplice_dna = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['raw_pred']['absplice']['absplice_preds']['dataset_splicemap']['dna'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        result_absplice_dna_variant = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['pred_variant_level']['absplice']['dataset_splicemap']['dna'],
    script:
        "../variant_level/variant_absplice_dna.py"
        
        
rule splicing_result_variant_level_absplice_rna_dataset_splicemap:
    input:
        pred_absplice_rna = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['raw_pred']['absplice']['absplice_preds']['dataset_splicemap']['rna'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        result_absplice_rna_variant = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['pred_variant_level']['absplice']['dataset_splicemap']['rna'],
    script:
        "../variant_level/variant_absplice_rna.py"