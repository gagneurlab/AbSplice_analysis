OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']

import pandas as pd
##==================================infer cats========================================= 
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