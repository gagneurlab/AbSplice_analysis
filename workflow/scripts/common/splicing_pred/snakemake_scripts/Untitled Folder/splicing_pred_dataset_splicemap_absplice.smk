OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
ABSPLICE_INPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['combined_input']

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