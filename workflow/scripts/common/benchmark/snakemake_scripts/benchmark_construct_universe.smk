##==================================BENCHMARK=========================================
rule benchmark_construct_universe:
    input:
        rare_variants = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['with_dist_to_exons'], #already filtered for protein coding
        count_table = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['updated'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1,
    output:
        rare_vars_tissue = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['rare_vars_tissue'],
        universe = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['universe'],
    script:
        "../construct_universe.py"
        
        
rule benchmark_construct_universe_variant_level:
    input:
        universe = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['universe'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1,
    output:
        universe = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['universe_variant_level'],
    script:
        "../construct_universe_variant_level.py"