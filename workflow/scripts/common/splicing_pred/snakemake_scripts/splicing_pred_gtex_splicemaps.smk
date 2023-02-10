OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
ABSPLICE_INPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['combined_input']

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
        
        
rule splicing_pred_infer_cat_gtex_splicemaps:
    input:
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
        mmsplice_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap']['gtex_splicemaps'],
        splicemap_5 = expand(config_precomputed['splicemap_gtex']['psi5'],
                             gtex_version=config['gtex_version'], tissue=config['gtex_tissues']),
        splicemap_3 = expand(config_precomputed['splicemap_gtex']['psi3'],
                             gtex_version=config['gtex_version'], tissue=config['gtex_tissues']),
        cat_count_table = expand(OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['updated'],
                                 tissue='{tissue_cat}')[0],
    params:
        tissue_cat = '{tissue_cat}',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        result = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['single_cat']
    script:
        "../infer_cat.py"
        

rule splicing_pred_mmsplice_splicemap_cat_combine:
    input:
        pred_mmsplice_cat = expand(OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['single_cat'],
                                   tissue_cat=config['tissues_cat'], vcf_id='{vcf_id}'),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        threads = 1
    output:
        result = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['all_cats']
    run:
        df_mmsplice_cat = pd.concat([pd.read_parquet(i).reset_index() for i in input.pred_mmsplice_cat])
        df_mmsplice_cat.to_parquet(output.result, index=False)
        
        
        
# --------------------------------AbSplice input------------------------------- 
rule splicing_pred_absplice_dna_input_gtex_splicemaps:
    input:
        pred_mmsplice = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap']['gtex_splicemaps'],
        pred_spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
        gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
        gene_tpm = config_precomputed['gene_tpm_gtex'].format(gtex_version=config['gtex_version']),
        var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
    script:
        "../absplice_dna_input.py"
        
        
rule splicing_pred_absplice_rna_input_gtex_splicemaps:
    input:
        absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
        pred_mmsplice_cat = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['single_cat'],
        CAT_pval = expand(OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['tissue_cat_pval'],
                          tissue='{tissue_cat}', vcf_id='{vcf_id}'),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['single_cat'],
    script:
        "../absplice_rna_input.py"
         
        
rule splicing_pred_absplice_rna_input_gtex_splicemaps_all_cats:
    input:
        absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
        pred_mmsplice_cat = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['gtex_splicemaps']['all_cats'],
        CAT_pval = expand(OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['tissue_cat_pval'],
                          tissue=config['tissues_cat'], vcf_id='{vcf_id}'),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['all_cats'],
    script:
        "../absplice_rna_input.py"