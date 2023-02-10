import os
import yaml
import pandas as pd
from tqdm import tqdm
from absplice_scripts.data.gtex_v8 import tissue_map, tissue_map_gene_expr

ABSPLICE_INPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['combined_input']

#------------------------------------------- all CATs --------------------------------------
rule absplice_rna_training_input_all_cats:
    input:
        universe = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['universe'],
        outliers = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist_variant_level'],
        absplice_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['all_cats'],
    params:
        tissue = '{tissue}',
        vcf_id = '{vcf_id}',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        absplice_input = OUTPUT_DIR_SPLICING + config['absplice_training']['training_input']['rna']['all_cats']
    script:
        "../absplice_combine_training_input.py"
        
        
rule absplice_rna_training_input_subset_all_cats:
    input:
        absplice_input = OUTPUT_DIR_SPLICING + config['absplice_training']['training_input']['rna']['all_cats'],
    params:
        tissue = '{tissue}',
        tissues_cat = config['tissues_cat'],
        vcf_id = '{vcf_id}',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        absplice_input = OUTPUT_DIR_SPLICING + config['absplice_training']['training_input_subset']['rna']['all_cats']
    script:
        "./training_input_subset.py"
        
        
rule absplice_rna_input_complete_all_cats:
    input:
        absplice_input = expand(OUTPUT_DIR_SPLICING + config['absplice_training']['training_input_subset']['rna']['all_cats'],
                                tissue=config['gtex_tissues'], vcf_id=wildcard_vcf_id)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        absplice_input = OUTPUT_DIR_SPLICING + config['absplice_training']['training_input_complete']['rna']['all_cats']
    run:
        df = pd.concat([pd.read_csv(i) for i in tqdm(input.absplice_input)])
        df = df[df['tissue'] != df['tissue_cat']]
        df.to_csv(output.absplice_input, index=False)
    
    
rule absplice_rna_training_5fold_all_cats:
    input:
        absplice_input = OUTPUT_DIR_SPLICING + config['absplice_training']['training_input_complete']['rna']['all_cats'],
    params:
        classifier = '{classifier}',
        feature_string = '{feature_string}',
        abs_features = '{abs_features}',
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        gene_tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        absplice_5fold_crossval = OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['rna']['all_cats']['5_fold_crossval'],
    script:
        "../absplice_training_5fold.py"
        
        
rule absplice_rna_training_whole_GTEx_all_cats:
    input:
        absplice_input = OUTPUT_DIR_SPLICING + config['absplice_training']['training_input_complete']['rna']['all_cats'],
    params:
        classifier = '{classifier}',
        feature_string = '{feature_string}',
        abs_features = '{abs_features}',
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        gene_tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        absplice_whole_GTEx = OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['rna']['all_cats']['whole_GTEx'],
    script:
        "../absplice_training_whole_GTEx.py"
        
        
rule absplice_rna_training_5fold_gene_level_all_cats:
    input:
        pred_absplice_rna = OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['rna']['all_cats']['5_fold_crossval'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_absplice_rna_gene = OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['rna']['all_cats']['gene_level'],
    script:
        "../gene_level_absplice_rna.py"
        
        
rule all_absplice_training_rna_all_cats:
    input:
        expand(rules.absplice_rna_training_5fold_gene_level_all_cats.output,
               classifier=config['absplice_params']['classifier'],
               abs_features=config['absplice_params']['abs_features'],
               feature_string=config['absplice_params']['feature_string_rna']
              ),
        expand(rules.absplice_rna_training_whole_GTEx_all_cats.output,
               classifier=config['absplice_params']['classifier'],
               abs_features=config['absplice_params']['abs_features'],
               feature_string=config['absplice_params']['feature_string_rna']
              ),