##==================================JUNCTION ANNOTATION=========================================
import pandas as pd

def analysis(wildcards):
    df_anno = pd.read_csv(config['DROP']['tissue_map'])
    tissue_map = dict(zip(df_anno['tissue'], df_anno['tissue_DROP']))
    return 'raw-' + tissue_map[wildcards['tissue']]

def fds(wildcards):
    return config['DROP']['splicemap']['working_dir'] + 'savedObjects/' + analysis(wildcards) + '/fds-object.RDS'


rule annotate_fraser:
    input:
        fds = fds
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000
    params:
        fraser_working_dir = config['DROP']['splicemap']['working_dir'],
        analysis = analysis
    output:
        raw_count_table = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['raw']
    script: 
        "../annotate_fraser.R"


rule annotate_count_table:
    input:
        count_table = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['raw'],
        gtf_file = config['gtf'],
        coding_genes = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['coding_genes'],
    params:
        tissue = '{tissue}',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
    output:
        count_table_with_annotation = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['with_annotation']
    script: 
        "../count_table_infer_annotation.py"
        
        
rule update_count_table:
    input:
        fasta = config['fasta'],
        raw_count_table = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['raw'],
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000
    params:
        tissue = '{tissue}',
        annotation_table = config['DROP']['sample_annotation'],
        update_samples = config['filtering_params']['count_table']['update_samples'],
        key_assay = config['filtering_params']['count_table']['key_assay'],
        value_assay = config['filtering_params']['count_table']['value_assay'],
#         remove_chr_from_chrom_annotation = config['filtering_params']['count_table']['remove_chr_from_chrom_annotation'],
        infer_strand = config['filtering_params']['count_table']['infer_strand'],
    output:
        updated_count_table = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['updated'],
    script:
        "../update_count_table.py"


rule all_count_table:
    input:
        expand(rules.annotate_fraser.output,
               tissue=config['splicemap_tissues']),
        expand(rules.annotate_count_table.output,
               tissue=config['splicemap_tissues']),
        expand(rules.update_count_table.output,
               tissue=config['splicemap_tissues']),