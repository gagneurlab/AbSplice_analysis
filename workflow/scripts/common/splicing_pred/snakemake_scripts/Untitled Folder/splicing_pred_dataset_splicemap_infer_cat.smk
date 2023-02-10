OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']

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