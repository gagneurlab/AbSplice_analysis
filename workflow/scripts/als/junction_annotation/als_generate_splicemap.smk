import pandas as pd
from absplice_scripts.data.als import tissue_map


rule als_update_sample_anno:
    input:
        sample_anno = config['DROP']['sample_annotation_pure'],
    output:
        sample_anno_updated = config['DROP']['sample_annotation'],
    script:
        "./correct_vcf_id_DROP_anno.py"
        
        
rule als_generate_gene_tpm:
    input:
        gene_count = config['gene_count'],
        gtf_file = config['gtf']
    params:
        tissue = 'motor_neuron',
    output:
        gene_count = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_count'],
        gene_tpm = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm'],
        gene_tpm_wide = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm_wide'],
    script:
        "./calculate_tpm.py"
        
        
rule als_generate_tissue_map:
    output:
        tissue_map = config['DROP']['tissue_map'],
    run:
        df = pd.DataFrame(tissue_map, index=[0]).transpose().reset_index().rename(columns={'index': 'tissue', 0: 'tissue_DROP'})
        df.to_csv(output.tissue_map, index=False)
        
        
rule als_generate_splicemap:
    input:
        count_table = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['updated'],
        gene_expression = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm_wide'],
        gtf_file = config['gtf']
    params:
        tissue = 'motor_neuron',
        method = config['method'],
        event_filter = config['event_filter'],
        percentile = config['filtering_params']['splicemap']['percentile'],
        percentile_min_read = config['filtering_params']['splicemap']['percentile_min_read'],
        median_cutoff = config['filtering_params']['splicemap']['median_cutoff'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 128000,
        ntasks = 1,
        threads = 10
    output:
        splicemap_psi5 = config_static['junction_annotation']['splicemap']['psi5'],
        splicemap_psi3 = config_static['junction_annotation']['splicemap']['psi3'],
    script:
        "./generate_splicemap.py"
        
        
rule all_als_generate_splicemap:
    input:
        expand(rules.als_generate_splicemap.output,
               tissue=config['splicemap_tissues'], genome=config['genome']),