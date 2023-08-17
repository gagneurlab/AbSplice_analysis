OUTPUT_DIR_SPLICING_POSTPROCESS = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['postprocess']
OUTPUT_DIR_SPLICING_AGG = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['variant_level']

groupby_index = ['variant', 'gene_id', 'sample']
groupby_index_tissue = ['variant', 'gene_id', 'sample', 'tissue']

# SQUIRLS
rule variant_level_squirls:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['squirls'],
    params:
        max_col = 'squirls_scores',
        groupby_index = groupby_index,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['squirls'],
    script:
        "./max_aggregation.py"

# SpliceAI
rule variant_level_spliceai:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['spliceai'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_score',
        groupby_index = groupby_index,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['spliceai']
    script:
        "./max_aggregation.py"

# SpliceAI + SpliceMap
rule variant_level_spliceai_splicemap:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['spliceai_splicemap']['all'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_score',
        groupby_index = groupby_index_tissue,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['spliceai_splicemap']['all']
    script:
        "./max_aggregation.py"

# SpliceAI + SpliceMap + PSI_ref
rule variant_level_spliceai_splicemap_ref_psi:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['all'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_score',
        groupby_index = groupby_index_tissue,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['all']
    script:
        "./max_aggregation.py"
        
# MMSplice
rule variant_level_mmsplice:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_logit_psi',
        groupby_index = groupby_index,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['mmsplice']
    script:
        "./max_aggregation.py"
        
# MTSplice      
rule variant_level_mtsplice:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mtsplice'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_logit_psi_mtsplice',
        groupby_index = groupby_index_tissue,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['mtsplice']
    script:
        "./max_aggregation.py"
        
# CADD-Splice      
rule variant_level_cadd_splice:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['cadd_splice'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'PHRED',
        groupby_index = groupby_index,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['cadd_splice']
    script:
        "./max_aggregation.py"

# MMSplice + SpliceMap
rule variant_level_mmsplice_splicemap:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice_splicemap'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_logit_psi',
        groupby_index = groupby_index_tissue,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['mmsplice_splicemap'],
    script:
        "./max_aggregation.py"

# MMSplice + SpliceMap  + PSI_ref  
rule variant_level_mmsplice_splicemap_ref_psi:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_psi',
        groupby_index = groupby_index_tissue,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi'],
    script:
        "./max_aggregation.py"

# CAT infer single CAT
rule variant_level_cat_infer_single_cat:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    params:
        max_col = 'delta_psi_cat',
        groupby_index = groupby_index_tissue,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat'],
    script:
        "./max_aggregation.py"

# CAT infer all CATs
rule variant_level_cat_infer_all_cats:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['all_cats'],
    params:
        max_col = 'delta_psi_cat',
        groupby_index = groupby_index_tissue,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['all_cats'],
    script:
        "./max_aggregation.py"

# AbSplice-DNA
rule variant_level_absplice_dna:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['absplice']['dna'],
    params:
        max_col = 'AbSplice_DNA',
        groupby_index = groupby_index_tissue,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['absplice']['dna'],
    script:
        "./max_aggregation.py"

# AbSplice-DNA (no samples)
rule variant_level_absplice_dna_no_samples:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['absplice']['dna_no_samples'],
    params:
        max_col = 'AbSplice_DNA',
        groupby_index = [x for x in groupby_index_tissue if x != 'sample'],
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['absplice']['dna_no_samples'],
    script:
        "./max_aggregation.py"

# AbSplice-DNA with extra info
rule variant_level_absplice_dna_extra_info:
    input:
        absplice = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['absplice']['dna'],
        spliceai = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['spliceai'],
        mmsplice_splicemap = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi'],
    output:
        absplice_extra_info = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['absplice_extra_info']['dna'],
    script:
        "./absplice_dna_add_extra_info.py"

# AbSplice-DNA with extra info
rule variant_level_absplice_dna_no_samples_extra_info:
    input:
        absplice = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['absplice']['dna_no_samples'],
        spliceai = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['spliceai'],
        mmsplice_splicemap = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi'],
    output:
        absplice_extra_info = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['absplice_extra_info']['dna_no_samples'],
    script:
        "./absplice_dna_add_extra_info.py"

# AbSplice-RNA single cat
rule variant_level_absplice_rna_single_cat:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['absplice']['rna']['single_cat'],
    params:
        max_col = 'AbSplice_RNA',
        groupby_index = groupby_index_tissue,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['absplice']['rna']['single_cat'],
    script:
        "./max_aggregation.py"

# AbSplice-RNA all cats
rule variant_level_absplice_rna_all_cats:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['absplice']['rna']['all_cats'],
    params:
        max_col = 'AbSplice_RNA',
        groupby_index = groupby_index_tissue,
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['absplice']['rna']['all_cats'],
    script:
        "./max_aggregation.py"



list_outputs = list()
if 'mmsplice' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_mmsplice.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'mtsplice' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_mtsplice.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'spliceai' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_spliceai.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'cadd_splice' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_cadd_splice.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'mmsplice_splicemap' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_mmsplice_splicemap.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'mmsplice_splicemap_ref_psi' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_mmsplice_splicemap_ref_psi.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'spliceai_splicemap' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_spliceai_splicemap.output, 
        vcf_id=wildcard_vcf_id), # SpliceAI + SpliceMap
    )
if 'spliceai_splicemap_ref_psi' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_spliceai_splicemap_ref_psi.output, 
        vcf_id=wildcard_vcf_id), # SpliceAI + SpliceMap + PSI_ref
    )
if 'cat_infer_single' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_cat_infer_single_cat.output, 
        vcf_id=wildcard_vcf_id, tissue_cat=config['tissues_cat']), # CAT infer (single CAT)
    )
if 'cat_infer_all' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_cat_infer_all_cats.output, 
        vcf_id=wildcard_vcf_id), # CAT infer (all CAT)
    )
if 'absplice_dna' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_absplice_dna.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-DNA
    )
if 'absplice_dna_no_samples' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_absplice_dna_no_samples.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-DNA (input)
    )
if 'absplice_rna_single' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_absplice_rna_single_cat.output, 
        vcf_id=wildcard_vcf_id, tissue_cat=config['tissues_cat']), # AbSplice-RNA single CAT
    )
if 'absplice_rna_all' in config['models']:
    list_outputs.append(
        expand(rules.variant_level_absplice_rna_all_cats.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-RNA all CATs
    )


rule all_splicing_result_variant_level:
    input:
        list_outputs


del OUTPUT_DIR_SPLICING_AGG
del groupby_index
del groupby_index_tissue