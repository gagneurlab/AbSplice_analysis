##==================================CODING GENES=========================================
rule gtf_junctions:
    input:
        gtf = config['gtf']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        threads = 1,
    output:
        gtf_junctions = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gtf_junctions']
    script:
        "../gtf_junctions.py"


rule coding_genes:
    input:
        gtf = config['gtf']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        threads = 1,
    output:
        coding_genes = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['coding_genes']
    script:
        "../coding_genes.py"

        
rule update_gene_mapping:
    input:
        gene_mapping_raw = config_precomputed['gene_id_to_name_mapping_raw'],
        coding_genes = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['coding_genes'],
    output:
        gene_mapping = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
    script:
        "../gene_mapping_fill_missing.py"