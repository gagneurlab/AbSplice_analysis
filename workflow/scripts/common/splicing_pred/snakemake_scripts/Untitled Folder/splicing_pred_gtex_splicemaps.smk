OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']

##==================================SPLICING PREDICTIONS=========================================     
rule splicing_pred_mmsplice_splicemap_gtex:
    input:
        vcf = config['vcf'],
        fasta = config['fasta'],
        splicemap_5 = expand(config_precomputed['splicemap_gtex']['psi5'],
                             gtex_version=config['gtex_version'], tissue=config['gtex_tissues']),
        splicemap_3 = expand(config_precomputed['splicemap_gtex']['psi3'],
                             gtex_version=config['gtex_version'], tissue=config['gtex_tissues']),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 4
    output:
        result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap']['gtex_splicemaps'])
    script:
        "../mmsplice_splicemap.py"
        

rule splicing_pred_spliceai_splicemap_gtex_tissue:
    input:
        spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
        splicemap_5 = expand(config_precomputed['splicemap_gtex']['psi5'],
                             gtex_version=config['gtex_version'], tissue='{tissue}'),
        splicemap_3 = expand(config_precomputed['splicemap_gtex']['psi3'],
                             gtex_version=config['gtex_version'], tissue='{tissue}'),
        gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
    params:
        tissue = '{tissue}'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 5000,
        threads = 1,
    output:
        spliceai_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['gtex_splicemaps']['tissue']
    script:
        "../spliceai_splicemap.py"

        
rule splicing_pred_spliceai_splicemap_gtex_all:
    input:
        df_tissue = expand(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['gtex_splicemaps']['tissue'],
                           vcf_id='{vcf_id}', tissue=config['gtex_tissues']),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1,
    output:
        df_all = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['gtex_splicemaps']['all']
    run:
        import pandas as pd
        from tqdm import tqdm
        df = pd.concat([pd.read_parquet(i) for i in tqdm(input.df_tissue)])
        df.to_parquet(output.df_all, index=False)
        
        
rule splicing_pred_spliceai_splicemap_ref_psi_gtex_tissue:
    input:
        spliceai_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['gtex_splicemaps']['tissue']
    params:
        tissue = '{tissue}',
        ref_psi_tolerance = 0.05
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1,
    output:
        spliceai_splicemap_ref_psi = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['gtex_splicemaps']['tissue']
    script:
        "../spliceai_splicemap_ref_psi.py"
        
        
rule splicing_pred_spliceai_splicemap_ref_psi_gtex_all:
    input:
        df_tissue = expand(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['gtex_splicemaps']['tissue'],
                           vcf_id='{vcf_id}', tissue=config['gtex_tissues']),
    output:
        df_all = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['gtex_splicemaps']['all']
    run:
        import pandas as pd
        from tqdm import tqdm
        df = pd.concat([pd.read_parquet(i) for i in tqdm(input.df_tissue)])
        df.to_parquet(output.df_all, index=False)