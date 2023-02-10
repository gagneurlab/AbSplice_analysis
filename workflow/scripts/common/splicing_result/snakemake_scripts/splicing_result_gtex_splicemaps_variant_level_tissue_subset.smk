OUTPUT_DIR_RAW = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
OUTPUT_DIR_TISSUE_SUBSET = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['tissue_subset'] + config_static['splicing_pred']['levels']['variant_level']

import pandas as pd            
# --------------------------------gene level tissue subset-------------------------------
rule splicing_result_tissue_subset_variant_level_spliceai_splicemap:
    input:
        pred_spliceai_splicemap = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['spliceai_splicemap']['gtex_splicemaps']['all'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
        gene_tpm = config_precomputed['gene_tpm_gtex'].format(gtex_version=config['gtex_version']),
        tissue_map = config_precomputed['gtex_tissue_map_main_tissue'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_spliceai_splicemap_gene = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['spliceai_splicemap']['gtex_splicemaps']['all']
    script:
        "../variant_level/gtex_tissue_subset/variant_spliceai_splicemap.py"
        
        
rule splicing_result_tissue_subset_variant_level_spliceai_splicemap_ref_psi:
    input:
        pred_spliceai_splicemap = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['gtex_splicemaps']['all'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
        gene_tpm = config_precomputed['gene_tpm_gtex'].format(gtex_version=config['gtex_version']),
        tissue_map = config_precomputed['gtex_tissue_map_main_tissue'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_spliceai_splicemap_gene = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['gtex_splicemaps']['all'],
    script:
        "../variant_level/gtex_tissue_subset/variant_spliceai_splicemap.py"
        
        
rule splicing_result_tissue_subset_variant_level_mmsplice:
    input:
        pred_mmsplice = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap']['gtex_splicemaps'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
        gene_tpm = config_precomputed['gene_tpm_gtex'].format(gtex_version=config['gtex_version']),
        tissue_map = config_precomputed['gtex_tissue_map_main_tissue'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_mmsplice_gene = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['mmsplice_splicemap']['gtex_splicemaps'],
    script:
        "../variant_level/gtex_tissue_subset/variant_mmsplice_splicemap.py"
        
        
rule splicing_result_tissue_subset_variant_level_mmsplice_ref_psi:
    input:
        pred_mmsplice = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap']['gtex_splicemaps'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
        gene_tpm = config_precomputed['gene_tpm_gtex'].format(gtex_version=config['gtex_version']),
        tissue_map = config_precomputed['gtex_tissue_map_main_tissue'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_mmsplice_gene = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi']['gtex_splicemaps'],
    script:
        "../variant_level/gtex_tissue_subset/variant_mmsplice_splicemap_ref_psi.py"
        
              
rule splicing_result_tissue_subset_variant_level_mmsplice_cat:
    input:
        pred_mmsplice_cat = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['single_cat'],
        tissue_map = config_precomputed['gtex_tissue_map_main_tissue'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_mmsplice_cat_gene = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['single_cat'],
    script:
        "../variant_level/gtex_tissue_subset/variant_mmsplice_splicemap_cat.py"
        
             
rule splicing_result_tissue_subset_variant_level_mmsplice_cat_all_cats:
    input:
        pred_mmsplice_cat = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['all_cats'],
        tissue_map = config_precomputed['gtex_tissue_map_main_tissue'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_mmsplice_cat_gene = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['all_cats'],
    script:
        "../variant_level/gtex_tissue_subset/variant_mmsplice_splicemap_cat.py"
              
        
rule splicing_result_tissue_subset_variant_level_absplice_dna:
    input:
        pred_absplice_dna = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
        tissue_map = config_precomputed['gtex_tissue_map_main_tissue'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_absplice_dna_gene = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
    script:
        "../variant_level/gtex_tissue_subset/variant_absplice_dna.py"
               
            
rule splicing_result_tissue_subset_variant_level_absplice_rna:
    input:
        pred_absplice_rna = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['single_cat'],
        tissue_map = config_precomputed['gtex_tissue_map_main_tissue'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_absplice_rna_gene = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['single_cat'],
    script:
        "../variant_level/gtex_tissue_subset/variant_absplice_rna.py"
        

rule splicing_result_tissue_subset_variant_level_absplice_rna_all_cats:
    input:
        pred_absplice_rna = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['all_cats'],
        tissue_map = config_precomputed['gtex_tissue_map_main_tissue'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_absplice_rna_gene = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['all_cats'],
    script:
        "../variant_level/gtex_tissue_subset/variant_absplice_rna.py"