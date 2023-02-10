def gnomad_version(wildcards):
    if wildcards['genome'] == 'hg38':
        return '3.1.2'
    elif wildcards['genome'] == 'hg19':
        return '2.1.1'
    else:
        raise ValueError('No such version. Use "hg38" or "hg19".')

        
rule download_gnomad_rocksdb:
    params:
        version = gnomad_version,
    resources:
        mem_mb = 8000,
        threads = 1
    output:
        output_path = directory(config_precomputed['gnomad']['maf_db']),
    shell:
        "gnomad_rocksdb_download --version {params.version} --db_path {output.output_path}"
