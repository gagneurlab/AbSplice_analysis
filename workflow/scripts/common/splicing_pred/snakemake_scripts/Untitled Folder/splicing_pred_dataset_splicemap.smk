OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']

##==================================SPLICING PREDICTIONS========================================= 
rule splicing_pred_mmsplice_splicemap_dataset:
    input:
        vcf = config['vcf'],
        fasta = config['fasta'],
        splicemap_5 = expand(OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['splicemap']['psi5'],
                            tissue=config['tissue_dataset']),
        splicemap_3 = expand(OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['splicemap']['psi3'],
                            tissue=config['tissue_dataset']),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 4
    output:
        result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap']['dataset_splicemap'])
    script:
        "../mmsplice_splicemap.py"
        
        
rule splicing_pred_spliceai_splicemap_dataset_tissue:
    input:
        spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
        splicemap_5 = expand(OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['splicemap']['psi5'],
                            tissue=config['tissue_dataset'], 
                             method=config['method'], event_filter=config['event_filter']),
        splicemap_3 = expand(OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['splicemap']['psi3'],
                            tissue=config['tissue_dataset'], 
                             method=config['method'], event_filter=config['event_filter']),
        gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
    params:
        tissue = '{tissue}'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1,
    output:
        spliceai_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['dataset_splicemap']['tissue']
    script:
        "../spliceai_splicemap.py"
        
        
rule splicing_pred_spliceai_splicemap_dataset_all:
    input:
        df_tissue = expand(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['dataset_splicemap']['tissue'],
                                    vcf_id='{vcf_id}', tissue=config['tissue_dataset']),
    output:
        df_all = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['dataset_splicemap']['all']
    run:
        import pandas as pd
        from tqdm import tqdm
        df = pd.concat([pd.read_parquet(i) for i in tqdm(input.df_tissue)])
        df.to_parquet(output.df_all, index=False)
        
        
rule splicing_pred_spliceai_splicemap_ref_psi_dataset_tissue:
    input:
        spliceai_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['dataset_splicemap']['tissue']
    params:
        tissue = '{tissue}',
        ref_psi_tolerance = 0.05
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1,
    output:
        spliceai_splicemap_ref_psi = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['dataset_splicemap']['tissue']
    script:
        "../spliceai_splicemap_ref_psi.py"
        
        
rule splicing_pred_spliceai_splicemap_ref_psi_dataset_all:
    input:
        df_tissue = expand(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['dataset_splicemap']['tissue'],
                                    vcf_id='{vcf_id}', tissue=config['tissue_dataset']),
    output:
        df_all = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['dataset_splicemap']['all']
    run:
        import pandas as pd
        from tqdm import tqdm
        df = pd.concat([pd.read_parquet(i) for i in tqdm(input.df_tissue)])
        df.to_parquet(output.df_all, index=False)