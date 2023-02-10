from absplice_scripts.utils.model_utils import *

OUTPUT_DIR_TISSUE_SUBSET = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['tissue_subset'] + config_static['splicing_pred']['levels']['gene_level']
OUTPUT_DIR_GENE = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['gene_level']

# -----------------------------------GTEx SpliceMaps-------------------------------------
rule benchmark_combine_dna_gtex_splicemaps:
    input:
        universe = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['universe'],
        outliers = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist_gene_level'],
        
        spliceai = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['spliceai'],
        spliceai_splicemap = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['spliceai_splicemap']['gtex_splicemaps']['all'],
        spliceai_splicemap_ref_psi = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['gtex_splicemaps']['all'],
        mmsplice_splicemap = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['mmsplice_splicemap']['gtex_splicemaps'],
        mmsplice_splicemap_ref_psi = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi']['gtex_splicemaps'],
        absplice_dna = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
        mmsplice = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['mmsplice'],
        mtsplice = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['mtsplice'],
        cadd_splice = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['cadd_splice'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    params:
        pred_tools = pred_tools_dna,
        tissue_specific_tools = tissue_specific_tools,
        splicemaps = 'gtex',
        gtex_tissue = '{gtex_tissue}',
        tissue = '{tissue}',
        cols_spliceai = [
            'variant', 'delta_score'],
        cols_mmsplice = [
            'variant', 'delta_logit_psi'],
        cols_mmsplice_splicemap = [
            'variant', 'delta_logit_psi', 'delta_psi', 'median_n', 'ref_psi', 'gene_tpm', 'tissue'],
        cols_mmsplice_splicemap_ref_psi = [
            'variant', 'delta_logit_psi', 'delta_psi', 'median_n', 'ref_psi', 'gene_tpm', 'tissue'],
        cols_absplice_dna = [
            'variant', 'AbSplice_DNA', 'tissue', 
            'delta_logit_psi', 'delta_psi', 'delta_score', 'splice_site_is_expressed'],
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
        combined_benchmark = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['combined_benchmark']['gtex_splicemaps']['dna'],
    script:
        "../benchmark_combine_dna.py"
        
        
rule performance_dna_gtex_splicemaps_tissue_subset:
    input:
        benchmark = expand(OUTPUT_DIR_BENCHMARK + config_static['benchmark']['combined_benchmark']['gtex_splicemaps']['dna'],
                           vcf_id=wildcard_vcf_id, tissue='{tissue}', gtex_tissue='{gtex_tissue}'),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    params:
        outlier_column = 'outlier',
        model_dict = model_dict_dna,
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        gene_tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    output:
        df_performance = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['performance']['gtex_splicemaps']['dna']['df'],
        aps_performance = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['performance']['gtex_splicemaps']['dna']['aps']
    script:
        "../performance_dna.py"
        
        
rule all_benchmark_dna_gtex_splicemaps:
    input:
        expand(rules.performance_dna_gtex_splicemaps_tissue_subset.output,
               tissue=config['tissue_dataset'],
               gtex_tissue=config['gtex_tissues'] 
              ),
        
        
