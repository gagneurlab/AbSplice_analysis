OUTPUT_DIR_SPLICING_POSTPROCESS = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['postprocess']
OUTPUT_DIR_SPLICING_AGG = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['gene_level']

groupby_index = ['gene_id', 'sample']
groupby_index_tissue = ['gene_id', 'sample', 'tissue']

# SQUIRLS
rule gene_level_squirls:
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
rule gene_level_spliceai:
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
rule gene_level_spliceai_splicemap:
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
rule gene_level_spliceai_splicemap_ref_psi:
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
rule gene_level_mmsplice:
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
rule gene_level_mtsplice:
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
rule gene_level_cadd_splice:
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
rule gene_level_mmsplice_splicemap:
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
rule gene_level_mmsplice_splicemap_ref_psi:
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
rule gene_level_cat_infer_single_cat:
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
rule gene_level_cat_infer_all_cats:
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
rule gene_level_absplice_dna:
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
rule gene_level_absplice_dna_no_samples:
    input:
        model = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['absplice']['dna_no_samples'],
    params:
        max_col = 'AbSplice_DNA',
        groupby_index = ['gene_id', 'tissue'],
    output:
        model_agg = OUTPUT_DIR_SPLICING_AGG + config_static['splicing_pred']['models']['absplice']['dna_no_samples'],
    script:
        "./max_aggregation.py"

# AbSplice-RNA single cat
rule gene_level_absplice_rna_single_cat:
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
rule gene_level_absplice_rna_all_cats:
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
        expand(rules.gene_level_mmsplice.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'mtsplice' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_mtsplice.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'spliceai' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_spliceai.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'cadd_splice' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_cadd_splice.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'mmsplice_splicemap' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_mmsplice_splicemap.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'mmsplice_splicemap_ref_psi' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_mmsplice_splicemap_ref_psi.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'spliceai_splicemap' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_spliceai_splicemap.output, 
        vcf_id=wildcard_vcf_id), # SpliceAI + SpliceMap
    )
if 'spliceai_splicemap_ref_psi' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_spliceai_splicemap_ref_psi.output, 
        vcf_id=wildcard_vcf_id), # SpliceAI + SpliceMap + PSI_ref
    )
if 'cat_infer_single' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_cat_infer_single_cat.output, 
        vcf_id=wildcard_vcf_id, tissue_cat=config['tissues_cat']), # CAT infer (single CAT)
    )
if 'cat_infer_all' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_cat_infer_all_cats.output, 
        vcf_id=wildcard_vcf_id), # CAT infer (all CAT)
    )
if 'absplice_dna' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_absplice_dna.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-DNA
    )
if 'absplice_dna_no_samples' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_absplice_dna_no_samples.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-DNA (input)
    )
if 'absplice_rna_single' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_absplice_rna_single_cat.output, 
        vcf_id=wildcard_vcf_id, tissue_cat=config['tissues_cat']), # AbSplice-RNA single CAT
    )
if 'absplice_rna_all' in config['models']:
    list_outputs.append(
        expand(rules.gene_level_absplice_rna_all_cats.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-RNA all CATs
    )


rule all_splicing_result_gene_level:
    input:
        list_outputs


del OUTPUT_DIR_SPLICING_AGG
del groupby_index
del groupby_index_tissue