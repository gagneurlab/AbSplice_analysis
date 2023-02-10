OUTPUT_DIR_SPLICING_RAW = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
OUTPUT_DIR_SPLICING_GENE_LEVEL = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['gene_level']
OUTPUT_DIR_SPLICING_VAR_LEVEL = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['variant_level']

##==================================SPLICING RESULTS=========================================
# rule squirls_anno:
#     input:
#         gtf = config['gtf'],
#         squirls = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['models']['squirls']['raw'],
#     output:
#         squirls_anno = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['models']['squirls']['with_gene_anno'],
#     script:
#         "../squirls_gene_anno.py"

# -------------------------------------gene level-------------------------------------
rule splicing_result_gene_level_spliceai:
    input:
        pred_spliceai = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['spliceai'],
        gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        result_spliceai_gene = OUTPUT_DIR_SPLICING_GENE_LEVEL + config_static['splicing_pred']['models']['spliceai']
    script:
        "../gene_level/gene_spliceai.py"
        

rule splicing_result_gene_level_mmsplice_baseline:
    input:
        pred_mmsplice_baseline = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['mmsplice'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        result_mmsplice_baseline_gene = OUTPUT_DIR_SPLICING_GENE_LEVEL + config_static['splicing_pred']['models']['mmsplice']
    script:
        "../gene_level/gene_mmsplice_baseline.py"
        
        
rule splicing_result_gene_level_mtsplice:
    input:
        pred_mtsplice = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['mtsplice'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        result_mtsplice_gene = OUTPUT_DIR_SPLICING_GENE_LEVEL + config_static['splicing_pred']['models']['mtsplice']
    script:
        "../gene_level/gene_mtsplice.py"
        
        
rule splicing_result_gene_level_cadd_splice:
    input:
        pred_cadd_splice = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['cadd_splice'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'PHRED'
    output:
        result_cadd_splice_gene = OUTPUT_DIR_SPLICING_GENE_LEVEL + config_static['splicing_pred']['models']['cadd_splice']
    script:
        "../gene_level/gene_cadd_splice.py"
        
        
# rule splicing_result_gene_level_squirls:
#     input:
#         pred_squirls = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['squirls']['with_gene_anno'],
#         var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 64000,
#         threads = 1
#     output:
#         result_squirls_gene = OUTPUT_DIR_SPLICING_GENE_LEVEL + config_static['splicing_pred']['models']['squirls']
#     script:
#         "../gene_level/gene_squirls.py"
        
        
# -------------------------------------variant level-------------------------------------
rule splicing_result_variant_level_spliceai:
    input:
        pred_spliceai = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['spliceai'],
        gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        result_spliceai_variant = OUTPUT_DIR_SPLICING_VAR_LEVEL + config_static['splicing_pred']['models']['spliceai']
    script:
        "../variant_level/variant_spliceai.py"
        

rule splicing_result_variant_level_mmsplice_baseline:
    input:
        pred_mmsplice_baseline = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['mmsplice'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        result_mmsplice_baseline_variant = OUTPUT_DIR_SPLICING_VAR_LEVEL + config_static['splicing_pred']['models']['mmsplice']
    script:
        "../variant_level/variant_mmsplice_baseline.py"
        
        
rule splicing_result_variant_level_mtsplice:
    input:
        pred_mtsplice = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['mtsplice'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        result_mtsplice_variant = OUTPUT_DIR_SPLICING_VAR_LEVEL + config_static['splicing_pred']['models']['mtsplice']
    script:
        "../variant_level/variant_mtsplice.py"
        
        
rule splicing_result_variant_level_cadd_splice:
    input:
        pred_cadd_splice = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['cadd_splice'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'PHRED'
    output:
        result_cadd_splice_variant = OUTPUT_DIR_SPLICING_VAR_LEVEL + config_static['splicing_pred']['models']['cadd_splice']
    script:
        "../variant_level/variant_cadd_splice.py"
        
        
# rule splicing_result_variant_level_squirls:
#     input:
#         pred_squirls = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['squirls']['with_gene_anno'],
#         var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 64000,
#         threads = 1
#     output:
#         result_squirls_variant = OUTPUT_DIR_SPLICING_VAR_LEVEL + config_static['splicing_pred']['models']['squirls']
#     script:
#         "../variant_level/variant_squirls.py"