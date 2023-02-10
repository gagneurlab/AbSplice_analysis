##==================================VCF ANNOTATION=========================================
rule vcf_annotate:
    input:
        vcf = config['vcf'],
        maf = ancient(config_precomputed['gnomad']['maf_db'].format(genome=config['genome'])),
    params:
        format_fields = config['filtering_params']['vcf']['format_fields'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 4
    output:
        vcf_annotation = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['raw'],
    script:
        "../vcf_annotation.py"
        
        
rule qual_filter_variants:
    input:
        rare_variants = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['raw'],
    params:
        min_GQ = config['filtering_params']['vcf']['min_GQ'],
        min_DP_ALT = config['filtering_params']['vcf']['min_DP_ALT'],
        max_length_indel = config['filtering_params']['vcf']['filter_long_indels'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        rare_variants_filtered = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    script:
        "../variants_qual_filter.py"
        
        
rule annotate_variant_distance:
    input:
        variants = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
        gtf = config['gtf'],
    params:
        remove_chr = config['filtering_params']['vcf']['remove_chr'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        variants_annotated = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['with_dist_to_exons'],
    script:
        "../variants_gene_id_exon_dist.py"
    