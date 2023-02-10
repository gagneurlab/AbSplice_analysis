OUTPUT_DIR_RAW = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
OUTPUT_DIR_GENE = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['gene_level']

import pandas as pd            
# --------------------------------gene level-------------------------------
# SpliceAI SpliceMap
rule splicing_result_gene_level_spliceai_splicemap_gtex_splicemaps:
    input:
        pred_spliceai_splicemap = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['spliceai_splicemap']['gtex_splicemaps']['all'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_spliceai_splicemap_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['spliceai_splicemap']['gtex_splicemaps']['all']
    script:
        "../gene_level/gene_spliceai_splicemap.py"
        
        
rule splicing_result_gene_level_spliceai_splicemap_ref_psi_gtex_splicemaps:
    input:
        pred_spliceai_splicemap = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['gtex_splicemaps']['all'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_spliceai_splicemap_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['gtex_splicemaps']['all']
    script:
        "../gene_level/gene_spliceai_splicemap.py"


# MMSplice SpliceMap
rule splicing_result_gene_level_mmsplice_splicemap_gtex_splicemaps:
    input:
        pred_mmsplice = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap']['gtex_splicemaps'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
        gene_tpm = config_precomputed['gene_tpm_gtex'].format(gtex_version=config['gtex_version']),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_mmsplice_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['mmsplice_splicemap']['gtex_splicemaps']
    script:
        "../gene_level/gene_mmsplice_splicemap.py"
        

# MMSplice SpliceMap Psi_ref
rule splicing_result_gene_level_mmsplice_splicemap_ref_psi_gtex_splicemaps:
    input:
        pred_mmsplice = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap']['gtex_splicemaps'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
        gene_tpm = config_precomputed['gene_tpm_gtex'].format(gtex_version=config['gtex_version']),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_mmsplice_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi']['gtex_splicemaps']
    script:
        "../gene_level/gene_mmsplice_splicemap_ref_psi.py"
        
        
rule splicing_result_gene_level_mmsplice_cat_gtex_splicemaps:
    input:
        pred_mmsplice_cat = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['single_cat'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_mmsplice_cat_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['single_cat'],
    script:
        "../gene_level/gene_mmsplice_splicemap_cat.py"
        
        
rule splicing_result_gene_level_mmsplice_cat_gtex_splicemaps_all_cats:
    input:
        pred_mmsplice_cat = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['all_cats'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_mmsplice_cat_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['all_cats'],
    script:
        "../gene_level/gene_mmsplice_splicemap_cat.py"