##==================================SPLICING PREDICTIONS========================================= 

OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']

# MMSplice
rule splicing_pred_mmsplice_baseline:
    input:
        vcf = config['vcf'],
        fasta = config['fasta'],
        gtf = config['gtf'],
    params:
        genome = config['genome'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1
    output:
        result = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice']
    script:
        "../mmsplice_exon.py"
        

# MTSplice
rule splicing_pred_mtsplice:
    input:
        vcf = config['vcf'],
        fasta = config['fasta'],
        gtf = config['gtf'],
#     conda: 'env/mtsplice_absplice.yaml' #mtsplice_absplice
    params:
        genome = config['genome'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1
    output:
        result = OUTPUT_DIR + config_static['splicing_pred']['models']['mtsplice']
    script:
        "../mtsplice.py"
        

# SpliceAI
rule splicing_pred_spliceai:
    input:
        vcf = config['vcf'],
        db = config_precomputed['spliceai']['db'].format(
            genome=config['genome']),
        fasta = config['fasta'],
    params:
        lookup_only = config['spliceai']['lookup_only'],
        genome = config['genome'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1,
#         gpu = 1,
    output:
        result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'])
    script:
        "../spliceai.py"
        
        
# CADD-Splice        
genome_mapper = {
    'hg38': 'GRCh38',
    'hg19': 'GRCh37',
}

rule cadd_splice:
    input:
        vcf = config['vcf'],
    conda: '/opt/modules/i12g/anaconda/envs/cadd_splice'
    params:
        genome = genome_mapper[config['genome']],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1
    output:
        result = OUTPUT_DIR + config_static['splicing_pred']['models']['cadd_splice']
    shell:
        '''/data/nasif12/home_if12/wagnern/Projects/CADD-scripts/CADD.sh -a -g "{params.genome}" -o {output.result} {input.vcf}'''
        
        
        
rule all_splicing_pred_baseline_models:
    input:
        expand(rules.splicing_pred_mmsplice_baseline.output, vcf_id=wildcard_vcf_id),
        expand(rules.splicing_pred_spliceai.output, vcf_id=wildcard_vcf_id),
        expand(rules.cadd_splice.output, vcf_id=wildcard_vcf_id),