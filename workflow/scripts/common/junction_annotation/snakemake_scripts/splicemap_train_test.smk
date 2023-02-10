##==================================SPLICEMAP TRAIN TEST=========================================  
rule splicemap_train_test_junction_rare_variants_faster:
    input:
        rare_vars = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['raw'],
        count_table_updated = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['updated'],
    params:
        tissue = '{tissue}',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 128000,
        threads = 1,
    output:
        junctions_with_rare_variant = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['train_test']['rare'],
    script:
        "../junction_rare_variants.py"


rule splicemap_train_test_junction_private_variants_faster:
    input:
        rare_variants = expand(OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['train_test']['rare'],
                                   vcf_id=wildcard_vcf_id, tissue='{tissue}'),
    params:
        tissue = '{tissue}',
        max_num_samples = 2,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1,
    output:
        private_variants = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['train_test']['private'],
    script:
        "../junction_private_variants.py"
        
        
        
rule splicemap_train_test_junction_train_test_split_faster:
    input:
        count_table_updated = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['updated'],
        private_variants = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['train_test']['private'],
    params:
        tissue = '{tissue}',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1,
    output:
        train_test_split = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['train_test']['tissue']
    script:
        "../junction_train_test_split.py"
        
        
        

