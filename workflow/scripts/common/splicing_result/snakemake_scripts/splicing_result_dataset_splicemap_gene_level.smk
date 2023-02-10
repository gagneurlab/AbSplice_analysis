OUTPUT_DIR_RAW = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
OUTPUT_DIR_GENE = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['gene_level']

import pandas as pd            
# --------------------------------gene level-------------------------------
# SpliceAI SpliceMap
rule splicing_result_gene_level_spliceai_splicemap_dataset_splicemap:
    input:
        pred_spliceai_splicemap = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['spliceai_splicemap']['dataset_splicemap']['all'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_spliceai_splicemap_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['spliceai_splicemap']['dataset_splicemap']['all']
    script:
        "../gene_level/gene_spliceai_splicemap.py"
        
        
rule splicing_result_gene_level_spliceai_splicemap_ref_psi_dataset_splicemap:
    input:
        pred_spliceai_splicemap = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['dataset_splicemap']['all'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_spliceai_splicemap_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['dataset_splicemap']['all']
    script:
        "../gene_level/gene_spliceai_splicemap.py"
        
# MMSplice SpliceMap 
rule splicing_result_gene_level_mmsplice_splicemap_dataset_splicemap:
    input:
        pred_mmsplice = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap']['dataset_splicemap'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
        gene_tpm = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        result_mmsplice_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['mmsplice_splicemap']['dataset_splicemap']
    script:
        "../gene_level/gene_mmsplice_splicemap.py"
        
        
rule splicing_result_gene_level_mmsplice_splicemap_ref_psi_dataset_splicemap:
    input:
        pred_mmsplice = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap']['dataset_splicemap'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
        gene_tpm = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        result_mmsplice_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi']['dataset_splicemap']
    script:
        "../gene_level/gene_mmsplice_splicemap_ref_psi.py"
        
        
rule splicing_result_gene_level_mmsplice_cat_dataset_splicemap:
    input:
        pred_mmsplice_cat = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['dataset_splicemap']['single_cat'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        result_mmsplice_cat_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['dataset_splicemap']['single_cat'],
    script:
        "../gene_level/gene_mmsplice_cat.py"
        

rule splicing_result_gene_level_mmsplice_cat_dataset_splicemap_all_cats:
    input:
        pred_mmsplice_cat = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['dataset_splicemap']['all_cats'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        result_mmsplice_cat_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['dataset_splicemap']['all_cats'],
    script:
        "../gene_level/gene_mmsplice_cat.py"
        

rule splicing_result_gene_level_absplice_dna_dataset_splicemap:
    input:
        pred_absplice_dna = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['dna'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        result_absplice_dna_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['dna'],
    script:
        "../gene_level/gene_absplice_dna.py"
               
            
rule splicing_result_gene_level_absplice_rna_dataset_splicemap:
    input:
        pred_absplice_rna = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['rna']['single_cat'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        result_absplice_rna_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['rna']['single_cat'],
    script:
        "../gene_level/gene_absplice_rna.py"
        
        
rule splicing_result_gene_level_absplice_rna_dataset_splicemap_all_cats:
    input:
        pred_absplice_rna = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['rna']['all_cats'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        result_absplice_rna_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['rna']['all_cats'],
    script:
        "../gene_level/gene_absplice_rna.py"

            
    