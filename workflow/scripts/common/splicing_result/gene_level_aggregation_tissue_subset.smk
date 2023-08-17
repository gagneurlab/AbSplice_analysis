OUTPUT_DIR_SPLICING_POSTPROCESS = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['postprocess']
OUTPUT_DIR_SPLICING_AGG = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['tissue_subset'] + config_static['splicing_pred']['levels']['gene_level']

groupby_index_tissue = ['gene_id', 'sample', 'tissue']

# SpliceAI + SpliceMap
rule gene_level_spliceai_splicemap__tissue_subset:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['spliceai_splicemap']['all'],
        tissue_map = config['tissue_map_subset'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_score',
        groupby_index = groupby_index_tissue,
        tissue_subset = True,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['spliceai_splicemap']['all']
    script:
        "./max_aggregation.py"

# SpliceAI + SpliceMap + PSI_ref
rule gene_level_spliceai_splicemap_ref_psi__tissue_subset:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['all'],
        tissue_map = config['tissue_map_subset'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_score',
        groupby_index = groupby_index_tissue,
        tissue_subset = True,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['all']
    script:
        "./max_aggregation.py"
        
# MTSplice      
rule gene_level_mtsplice__tissue_subset:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mtsplice'],
        tissue_map = config['tissue_map_subset'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_logit_psi_mtsplice',
        groupby_index = groupby_index_tissue,
        tissue_subset = True,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['mtsplice']
    script:
        "./max_aggregation.py"
        
# MMSplice + SpliceMap
rule gene_level_mmsplice_splicemap__tissue_subset:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice_splicemap'],
        tissue_map = config['tissue_map_subset'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_logit_psi',
        groupby_index = groupby_index_tissue,
        tissue_subset = True,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['mmsplice_splicemap'],
    script:
        "./max_aggregation.py"

# MMSplice + SpliceMap  + PSI_ref  
rule gene_level_mmsplice_splicemap_ref_psi__tissue_subset:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi'],
        tissue_map = config['tissue_map_subset'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_psi',
        groupby_index = groupby_index_tissue,
        tissue_subset = True,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi'],
    script:
        "./max_aggregation.py"

# CAT infer single CAT
rule gene_level_cat_infer_single_cat__tissue_subset:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat'],
        tissue_map = config['tissue_map_subset'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_psi_cat',
        groupby_index = groupby_index_tissue,
        tissue_subset = True,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat'],
    script:
        "./max_aggregation.py"

# CAT infer all CATs
rule gene_level_cat_infer_all_cats__tissue_subset:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['all_cats'],
        tissue_map = config['tissue_map_subset'],
    params:
        max_col = 'delta_psi_cat',
        groupby_index = groupby_index_tissue,
        tissue_subset = True,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['all_cats'],
    script:
        "./max_aggregation.py"

# AbSplice-DNA
rule gene_level_absplice_dna__tissue_subset:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['absplice']['dna'],
        tissue_map = config['tissue_map_subset'],
    params:
        max_col = 'AbSplice_DNA',
        groupby_index = groupby_index_tissue,
        tissue_subset = True,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['absplice']['dna'],
    script:
        "./max_aggregation.py"

# AbSplice-DNA (no samples)
rule gene_level_absplice_dna_no_samples__tissue_subset:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['absplice']['dna_no_samples'],
        tissue_map = config['tissue_map_subset'],
    params:
        max_col = 'AbSplice_DNA',
        groupby_index = groupby_index_tissue,
        tissue_subset = True,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['absplice']['dna_no_samples'],
    script:
        "./max_aggregation.py"

# AbSplice-RNA single cat
rule gene_level_absplice_rna_single_cat__tissue_subset:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['absplice']['rna']['single_cat'],
        tissue_map = config['tissue_map_subset'],
    params:
        max_col = 'AbSplice_RNA',
        groupby_index = groupby_index_tissue,
        tissue_subset = True,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['absplice']['rna']['single_cat'],
    script:
        "./max_aggregation.py"

# AbSplice-RNA all cats
rule gene_level_absplice_rna_all_cats__tissue_subset:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['absplice']['rna']['all_cats'],
        tissue_map = config['tissue_map_subset'],
    params:
        max_col = 'AbSplice_RNA',
        groupby_index = groupby_index_tissue,
        tissue_subset = True,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['absplice']['rna']['all_cats'],
    script:
        "./max_aggregation.py"


list_outputs = list()
if 'mtsplice' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_mtsplice__tissue_subset.output, 
        vcf_id=wildcard_vcf_id, tissue_pred=config['tissues_subset']), # MMSplice + SpliceMap
    )
if 'mmsplice_splicemap' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_mmsplice_splicemap__tissue_subset.output, 
        vcf_id=wildcard_vcf_id, tissue_pred=config['tissues_subset']), # MMSplice + SpliceMap
    )
if 'mmsplice_splicemap_ref_psi' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_mmsplice_splicemap_ref_psi__tissue_subset.output, 
        vcf_id=wildcard_vcf_id, tissue_pred=config['tissues_subset']), # MMSplice + SpliceMap
    )
if 'spliceai_splicemap' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_spliceai_splicemap__tissue_subset.output, 
        vcf_id=wildcard_vcf_id, tissue_pred=config['tissues_subset']), # SpliceAI + SpliceMap
    )
if 'spliceai_splicemap_ref_psi' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_spliceai_splicemap_ref_psi__tissue_subset.output, 
        vcf_id=wildcard_vcf_id, tissue_pred=config['tissues_subset']), # SpliceAI + SpliceMap + PSI_ref
    )
if 'cat_infer_single' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_cat_infer_single_cat__tissue_subset.output, 
        vcf_id=wildcard_vcf_id, tissue_pred=config['tissues_subset'], tissue_cat=config['tissues_cat']), # CAT infer (single CAT)
    )
if 'cat_infer_all' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_cat_infer_all_cats__tissue_subset.output, 
        vcf_id=wildcard_vcf_id, tissue_pred=config['tissues_subset']), # CAT infer (all CAT)
    )
if 'absplice_dna' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_absplice_dna__tissue_subset.output, 
        vcf_id=wildcard_vcf_id, tissue_pred=config['tissues_subset']), # AbSplice-DNA
    )
if 'absplice_dna_no_samples' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_absplice_dna_no_samples__tissue_subset.output, 
        vcf_id=wildcard_vcf_id, tissue_pred=config['tissues_subset']), # AbSplice-DNA (input)
    )
if 'absplice_rna_single' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_absplice_rna_single_cat__tissue_subset.output, 
        vcf_id=wildcard_vcf_id, tissue_pred=config['tissues_subset'], tissue_cat=config['tissues_cat']), # AbSplice-RNA single CAT
    )
if 'absplice_rna_all' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_absplice_rna_all_cats__tissue_subset.output, 
        vcf_id=wildcard_vcf_id, tissue_pred=config['tissues_subset']), # AbSplice-RNA all CATs
    )


rule all_splicing_result_gene_level__tissue_subset:
    input:
        list_outputs

del OUTPUT_DIR_SPLICING_AGG
del groupby_index_tissue