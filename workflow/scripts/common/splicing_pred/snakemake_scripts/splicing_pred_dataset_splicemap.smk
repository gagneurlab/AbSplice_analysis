OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
ABSPLICE_INPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['combined_input']

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
        
        
        
##==================================infer cats=========================================  
rule splicing_result_infer_cat_dataset_splicemap:
    input:
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
        mmsplice_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap']['dataset_splicemap'],
        splicemap_5 = expand(OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['splicemap']['psi5'],
                            tissue='{tissue}', method='{method}', event_filter='{event_filter}'),
        splicemap_3 = expand(OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['splicemap']['psi3'],
                            tissue='{tissue}', method='{method}', event_filter='{event_filter}'),
        cat_count_table = expand(OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['updated'],
                                 tissue='{tissue_cat}'),
    params:
        tissue_cat = '{tissue_cat}',
    resources:
        mem_mb = 128000,
        threads = 4
    output:
        result = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['dataset_splicemap']['single_cat']
    script:
        "../infer_cat.py"
        

##==================================AbSplice=========================================  
rule splicing_pred_absplice_dna_input_dataset_splicemap:
    input:
        pred_mmsplice = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap']['dataset_splicemap'],
        pred_spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
        gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
        gene_tpm = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm'],
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['dna'],
    script:
        "../absplice_dna_input.py"
         
        
rule splicing_pred_absplice_rna_input_dataset_splicemap:
    input:
        absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['dna'],
        pred_mmsplice_cat = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['dataset_splicemap']['single_cat'],
        # cat_pval = TODO
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['rna']['single_cat'],
    script:
        "../absplice_rna_input.py"
        
        
rule splicing_pred_absplice_rna_input_dataset_splicemap_all_cats:
    input:
        absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['dna'],
        pred_mmsplice_cat = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['dataset_splicemap']['all_cats'],
        # cat_pval = TODO
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['rna']['all_cats'],
    script:
        "../absplice_rna_input.py"
        
        
# --------------------------------predict AbSplice------------------------------- 
rule splicing_pred_absplice_dna_dataset_splicemap:
    input:
        absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['dna'],
        absplice_model = config_precomputed['absplice']['dna'],
    params:
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        absplice_dna_pred = OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['dna'],
    script:
        "../absplice_dna.py"
        
        
rule splicing_pred_absplice_rna_dataset_splicemap:
    input:
        absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['rna']['single_cat'],
        absplice_model = config_precomputed['absplice']['rna'],
    params:
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        absplice_rna_pred = OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['rna']['single_cat'],
    script:
        "../absplice_rna.py"
        
        
rule splicing_pred_absplice_rna_dataset_splicemap_all_cats:
    input:
        absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['rna']['all_cats'],
        absplice_model = config_precomputed['absplice']['rna'],
    params:
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    resources:
        mem_mb = 32000,
        threads = 4
    output:
        absplice_rna_pred = OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['dataset_splicemap']['rna']['all_cats'],
    script:
        "../absplice_rna.py"