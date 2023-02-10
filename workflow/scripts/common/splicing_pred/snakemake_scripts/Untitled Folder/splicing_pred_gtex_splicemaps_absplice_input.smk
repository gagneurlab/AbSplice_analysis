OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
ABSPLICE_INPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['combined_input']

import pandas as pd
# --------------------------------AbSplice input------------------------------- 
rule splicing_pred_absplice_dna_input_gtex_splicemaps:
    input:
        pred_mmsplice = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap']['gtex_splicemaps'],
        pred_spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
        gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
        gene_tpm = config_precomputed['gene_tpm_gtex'].format(gtex_version=config['gtex_version']),
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
    script:
        "../absplice_dna_input.py"
        
        
rule splicing_pred_absplice_rna_input_gtex_splicemaps:
    input:
        absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
        pred_mmsplice_cat = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['single_cat'],
        CAT_pval = expand(OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['tissue_cat_pval'],
                          tissue='{tissue_cat}', vcf_id='{vcf_id}'),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['single_cat'],
    script:
        "../absplice_rna_input.py"
         
        
rule splicing_pred_absplice_rna_input_gtex_splicemaps_all_cats:
    input:
        absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
        pred_mmsplice_cat = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['all_cats'],
        CAT_pval = expand(OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['tissue_cat_pval'],
                          tissue=config['tissues_cat'], vcf_id='{vcf_id}'),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['all_cats'],
    script:
        "../absplice_rna_input.py"