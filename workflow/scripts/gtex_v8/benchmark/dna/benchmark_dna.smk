from absplice_scripts.utils.model_utils import *

OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['gene_level']

# -----------------------------------GTEx SpliceMaps-------------------------------------
rule benchmark_combine_dna_gtex_splicemaps_GTEx:
    input:
        universe = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['universe'],
        outliers = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist_gene_level'],
        
        spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
        spliceai_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['gtex_splicemaps']['all'],
        spliceai_splicemap_ref_psi = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['gtex_splicemaps']['all'],
        mmsplice_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap']['gtex_splicemaps'],
        mmsplice_splicemap_ref_psi = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi']['gtex_splicemaps'],
        absplice_dna = expand(OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['dna']['gene_level'], 
                              classifier='{classifier}', feature_string='{feature_string_dna}', abs_features='{abs_features}')[0],
        mmsplice = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice'],
        mtsplice = OUTPUT_DIR + config_static['splicing_pred']['models']['mtsplice'],
        cadd_splice = OUTPUT_DIR + config_static['splicing_pred']['models']['cadd_splice'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    params:
        pred_tools = pred_tools_dna,
        tissue_specific_tools = tissue_specific_tools,
        unique_index = ['gene_id', 'sample'],
        cols_spliceai = [
            'variant', 'delta_score'],
        cols_mmsplice = [
            'variant', 'delta_logit_psi'],
        cols_mmsplice_splicemap = [
            'variant', 'delta_logit_psi', 'delta_psi', 'median_n', 'ref_psi', 'gene_tpm', 'tissue'],
        cols_mmsplice_splicemap_ref_psi = [
            'variant', 'delta_logit_psi', 'delta_psi', 'median_n', 'ref_psi', 'gene_tpm', 'tissue'],
        cols_absplice = [
            'variant', 'AbSplice_DNA', 'tissue'],
        cols_mtsplice = [
            'variant', 'delta_logit_psi_mtsplice', 'tissue'
        ],
        cols_cadd_splice = [
            'variant', 'PHRED'
        ],
        cols_spliceai_splicemap = [
            'junctions', 'delta_score', 'psi5_median_n', 'psi3_median_n', 'psi5_ref_psi', 'psi3_ref_psi', 'tissue'
        ],
        cols_spliceai_splicemap_ref_psi = [
            'junctions', 'delta_score', 'psi5_median_n', 'psi3_median_n', 'psi5_ref_psi', 'psi3_ref_psi', 'tissue'
        ]
    output:
        combined_benchmark = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['combined_benchmark']['gtex_splicemaps']['dna'],
    script:
        "../benchmark_combine.py"
        
            
rule performance_dna_all_tissues:
    input:
        benchmark = expand(OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['combined_benchmark']['gtex_splicemaps']['dna'],
                           vcf_id=wildcard_vcf_id, tissue=config['gtex_tissues'],
                           classifier='{classifier}', feature_string_dna='{feature_string_dna}', 
                           event_filter='{event_filter}', abs_features='{abs_features}'),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    params:
        subset_tissue = False,
        outlier_column = 'outlier',
        model_dict = model_dict_dna,
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        gene_tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    output:
        df_performance = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['performance']['gtex_splicemaps']['dna']['all_tissues']['df'],
        aps_performance = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['performance']['gtex_splicemaps']['dna']['all_tissues']['aps']
    script:
        "../../../../scripts/common/benchmark/performance_dna.py"
        
        
rule performance_dna_single_tissue:
    input:
        benchmark = expand(OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['combined_benchmark']['gtex_splicemaps']['dna'],
                           vcf_id=wildcard_vcf_id, tissue='{tissue}',
                           classifier='{classifier}', feature_string_dna='{feature_string_dna}', 
                           event_filter='{event_filter}', abs_features='{abs_features}'),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    params:
        # subset_tissue = True,
        outlier_column = 'outlier',
        model_dict = model_dict_dna,
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        gene_tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    output:
        df_performance = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['performance']['gtex_splicemaps']['dna']['single_tissue']['df'],
        aps_performance = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['performance']['gtex_splicemaps']['dna']['single_tissue']['aps']
    script:
        "../../../../scripts/common/benchmark/performance_dna.py"
        
        
rule performance_dna_across_tissues_boxplot:
    input:
        benchmark = expand(OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['performance']['gtex_splicemaps']['dna']['single_tissue']['aps'],
                           tissue=config['gtex_tissues'], classifier='{classifier}',
                           feature_string_dna='{feature_string_dna}', abs_features='{abs_features}'),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    output:
        df_performance_boxplot = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['performance']['gtex_splicemaps']['dna']['across_tissues']['boxplot_aps'],
    script:
        "../performance_across_tissues_boxplot.py"
        
        
rule corresponding_thresholds:
    input:
        aps_performance = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['performance']['gtex_splicemaps']['dna']['all_tissues']['aps']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    output:
        thresholds = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['performance']['gtex_splicemaps']['dna']['thresholds']['corresponding_thresholds'],
        thresholds_per_model = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['performance']['gtex_splicemaps']['dna']['thresholds']['thresholds_per_model'],
    script:
        "../corresponding_thresholds.py"
        
        
rule threshold_points_pr_curve:
    input:
        thresholds_per_model = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['performance']['gtex_splicemaps']['dna']['thresholds']['thresholds_per_model'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    output:
        threshold_points_pr_curve = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params'] + config['benchmark']['performance']['gtex_splicemaps']['dna']['thresholds']['threshold_points_pr_curve'],
    script:
        "../get_pr_for_thresholds.py"
        
        
rule all_benchmark_dna:
    input:
        expand(rules.performance_dna_all_tissues.output,
               classifier=config['absplice_params']['classifier'],
               abs_features=config['absplice_params']['abs_features'],
               feature_string_dna=config['absplice_params']['feature_string_dna'],
              ),
        expand(rules.performance_dna_across_tissues_boxplot.output,
               classifier=config['absplice_params']['classifier'],
               abs_features=config['absplice_params']['abs_features'],
               feature_string_dna=config['absplice_params']['feature_string_dna'],
              ),
        expand(rules.threshold_points_pr_curve.output,
               classifier=config['absplice_params']['classifier'],
               abs_features=config['absplice_params']['abs_features'],
               feature_string_dna=config['absplice_params']['feature_string_dna'],
              ),