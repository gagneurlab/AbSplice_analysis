import pandas as pd

include: "./junction_annotation_coding_genes.smk"
include: "./junction_annotation_count_table.smk"
include: "./splicemap_train_test.smk"

rule generate_splicemap:
    input:
        count_table_updated = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['updated'],
        gtf_file = config['gtf'],
        train_test_split = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['train_test']['tissue'],
        gene_expression = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm_wide'],
    params:
        tissue = '{tissue}',
        method = config['method'],
        event_filter = config['event_filter'],
        percentile = config['filtering_params']['splicemap']['percentile'],
        percentile_min_read = config['filtering_params']['splicemap']['percentile_min_read'],
        median_cutoff = config['filtering_params']['splicemap']['median_cutoff'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 128000,
        ntasks = 1,
        threads = 1
    output:
        splicemap_psi5 = config_static['junction_annotation']['splicemap']['psi5'],
        splicemap_psi3 = config_static['junction_annotation']['splicemap']['psi3'],
    script:
        "../generate_splicemap.py"
        
        
rule all_generate_splicemap:
    input:
        expand(rules.generate_splicemap.output,
               tissue=config['splicemap_tissues'], genome=config['genome']),
    