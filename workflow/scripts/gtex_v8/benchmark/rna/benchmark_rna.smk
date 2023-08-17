from absplice_scripts.utils.model_utils import *

OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['gene_level']

rule samples_with_cat:
    input:
        sample_anno = config['DROP']['sample_annotation'],
        tissue_map = config['DROP']['tissue_map'],
    params:
        tissues_cat = config['tissues_cat']
    output:
        samples_with_cat = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['samples_with_cat'],
    script:
        "./samples_with_cat.py"
        
        
rule universe_with_cat:
    input:
        universe = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['universe'],
        samples_with_cat = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['samples_with_cat'],
    output:
        universe_with_cat = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['universe_cat'],
    script:
        "./universe_with_cat.py"

        
# -----------------------------------GTEx SpliceMaps-------------------------------------
rule benchmark_combine_rna_gtex_splicemaps_GTEx:
    input:
        universe = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['universe_cat'],
        outliers = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist_gene_level'],
        
        spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
        spliceai_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['all'],
        spliceai_splicemap_ref_psi = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['all'],
        mmsplice_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap'],
        mmsplice_splicemap_ref_psi = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi'],
        absplice_dna = expand(OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['dna']['gene_level'], 
                              classifier='{classifier}', feature_string='{feature_string_dna}', abs_features='{abs_features}')[0],
        
        mmsplice_splicemap_cat = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat'],
        cat_pval = expand(OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['tissue_cat_pval'],
                          tissue='{tissue_cat}', vcf_id='{vcf_id}')[0],
        absplice_rna_single_cat = expand(OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['rna']['single_cat']['gene_level'],
                                         tissue_cat='{tissue_cat}',
                                         classifier='{classifier}', feature_string='{feature_string_rna}', abs_features='{abs_features}')[0],
        absplice_rna_all_cats = expand(OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['rna']['all_cats']['gene_level'],
                                       classifier='{classifier}', feature_string='{feature_string_rna}', abs_features='{abs_features}')[0],
        
        mmsplice = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice'],
        mtsplice = OUTPUT_DIR + config_static['splicing_pred']['models']['mtsplice'],
        cadd_splice = OUTPUT_DIR + config_static['splicing_pred']['models']['cadd_splice'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    params:
        pred_tools = pred_tools_rna,
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
        cols_absplice_dna = [
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
        ],
        
        
        cols_mmsplice_splicemap_cat = [
            'variant', 'delta_psi', 'median_n', 'ref_psi', 'tissue',
            'count_cat', 'median_n_cat', 'psi_cat', 'ref_psi_cat', 'delta_psi_cat', 'tissue_cat'],
        cols_absplice_rna_single_cat = [
            'AbSplice_RNA', 'tissue'],
        cols_absplice_rna_all_cats = [
            'AbSplice_RNA', 'tissue'],
        cols_cat_pval = [
            'pValueGene_g_minus_log10'],
    output:
        combined_benchmark = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['combined_benchmark']['rna'],
    script:
        "../benchmark_combine.py"
        
            
rule performance_rna_all_tissues:
    input:
        benchmark = expand(OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['combined_benchmark']['rna'],
                           vcf_id=wildcard_vcf_id, tissue=config['tissues'], tissue_cat='{tissue_cat}',
                           feature_string_dna='{feature_string_dna}', feature_string_rna='{feature_string_rna}', 
                           classifier='{classifier}', abs_features='{abs_features}'),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    params:
        subset_tissue = False,
        outlier_column = 'outlier',
        model_dict = model_dict_dna,
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        gene_tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    output:
        df_performance = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['performance']['rna']['all_tissues']['df'],
        aps_performance = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['performance']['rna']['all_tissues']['aps']
    script:
        "../../../../scripts/common/benchmark/performance_rna.py"
        
        
rule performance_rna_single_tissue:
    input:
        benchmark = expand(OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['combined_benchmark']['rna'],
                           vcf_id=wildcard_vcf_id, tissue='{tissue}', tissue_cat='{tissue_cat}',
                           feature_string_dna='{feature_string_dna}', feature_string_rna='{feature_string_rna}', 
                           classifier='{classifier}', abs_features='{abs_features}'),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    params:
        outlier_column = 'outlier',
        model_dict = model_dict_dna,
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        gene_tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    output:
        df_performance = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['performance']['rna']['single_tissue']['df'],
        aps_performance = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['performance']['rna']['single_tissue']['aps']
    script:
        "../../../../scripts/common/benchmark/performance_rna.py"
        
        
rule performance_rna_across_tissues_boxplot:
    input:
        benchmark = expand(OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['performance']['rna']['single_tissue']['aps'],
                           tissue=config['tissues'], tissue_cat='{tissue_cat}', cat_pairing='{cat_pairing}',
                           feature_string_dna='{feature_string_dna}', feature_string_rna='{feature_string_rna}', 
                           classifier='{classifier}', abs_features='{abs_features}'),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    output:
        df_performance_boxplot = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['performance']['rna']['across_tissues']['boxplot_aps'],
    script:
        "../performance_across_tissues_boxplot.py"
        
        
rule corresponding_thresholds_rna:
    input:
        aps_performance = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['performance']['rna']['all_tissues']['aps']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    output:
        thresholds = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['performance']['rna']['thresholds']['corresponding_thresholds'],
        thresholds_per_model = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['performance']['rna']['thresholds']['thresholds_per_model'],
    script:
        "../corresponding_thresholds.py"
        
        
rule threshold_points_pr_curve_rna:
    input:
        thresholds_per_model = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['performance']['rna']['thresholds']['thresholds_per_model'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    output:
        threshold_points_pr_curve = OUTPUT_DIR_BENCHMARK + config['benchmark']['absplice_model_params_rna'] + config['benchmark']['performance']['rna']['thresholds']['threshold_points_pr_curve'],
    script:
        "../get_pr_for_thresholds.py"


def filter_wildcards(path_list):
    correct_paths = list()
    for path in path_list:
        wildcard = re.findall('\/samples_paired\/(.+)\/tissue_cat=(.*?)\/', path)[0]
        if wildcard[1] in wildcard[0]:
            correct_paths.append(path)
    return correct_paths
    

rule all_benchmark_rna:
    input:
        filter_wildcards(
            expand(rules.performance_rna_across_tissues_boxplot.output,
                   classifier=config['absplice_params']['classifier'],
                   abs_features=config['absplice_params']['abs_features'],
                   feature_string_dna=config['absplice_params']['feature_string_dna'],
                   feature_string_rna=config['absplice_params']['feature_string_rna'],
                   cat_pairing=config['cat_pairing'],
                   tissue_cat=config['tissues_cat'],
                  ),
        ),
        filter_wildcards(
            expand(rules.performance_rna_all_tissues.output,
                   classifier=config['absplice_params']['classifier'],
                   abs_features=config['absplice_params']['abs_features'],
                   feature_string_dna=config['absplice_params']['feature_string_dna'],
                   feature_string_rna=config['absplice_params']['feature_string_rna'],
                   cat_pairing=config['cat_pairing'],
                   tissue_cat=config['tissues_cat'],
                  ),
        )