OUTPUT_DIR_SPLICING_RAW = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
OUTPUT_DIR_SPLICING_POSTPROCESS = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['postprocess']

var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered']

# MMSplice
rule splicing_result_postprocess_mmsplice:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['mmsplice'],
        var_samples_df = var_samples_df,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice']
    script:
        "./mmsplice_baseline.py"
        
# MTSplice      
rule splicing_result_postprocess_mtsplice:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['mtsplice'],
        var_samples_df = var_samples_df,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mtsplice']
    script:
        "./mtsplice.py"

# SQUIRLS
rule splicing_result_postprocess_squirls:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['squirls_raw'],
        gtf = config['gtf'],
        var_samples_df = var_samples_df,
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['squirls'],
    script:
        "./squirls.py"

# SpliceAI
rule splicing_result_postprocess_spliceai:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['spliceai'],
        gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
        var_samples_df = var_samples_df,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['spliceai']
    script:
        "./spliceai.py"

# SpliceAI + SpliceMap
rule splicing_result_postprocess_spliceai_splicemap:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['spliceai_splicemap']['all'],
        var_samples_df = var_samples_df,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['spliceai_splicemap']['all']
    script:
        "./spliceai_splicemap.py"

# SpliceAI + SpliceMap + PSI_ref
rule splicing_result_postprocess_spliceai_splicemap_ref_psi:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['all'],
        var_samples_df = var_samples_df,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['all']
    script:
        "./spliceai_splicemap.py"
        
# CADD-Splice      
rule splicing_result_postprocess_cadd_splice:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['cadd_splice'],
        var_samples_df = var_samples_df,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['cadd_splice']
    script:
        "./cadd_splice.py"

# MMSplice + SpliceMap
rule splicing_result_postprocess_mmsplice_splicemap:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap'],
        gene_tpm = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm'],
        var_samples_df = var_samples_df,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice_splicemap'],
    script:
        "./mmsplice_splicemap.py"

# MMSplice + SpliceMap  + PSI_ref  
rule splicing_result_postprocess_mmsplice_splicemap_ref_psi:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap'],
        gene_tpm = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm'],
        var_samples_df = var_samples_df,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi'],
    script:
        "./mmsplice_splicemap.py"

# ----------------- models that do not need postprocessing, just copy for convenience of structure ------------------
# CAT infer single CAT
rule splicing_result_postprocess_cat_infer_single_cat:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat'],
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat'],
    shell:
        'cp {input} {output}'

# CAT infer all CATs
rule splicing_result_postprocess_cat_infer_all_cats:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['all_cats'],
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['all_cats'],
    shell:
        'cp {input} {output}'

# AbSplice-DNA
rule splicing_result_postprocess_absplice_dna:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['absplice']['dna'],
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['absplice']['dna'],
    shell:
        'cp {input} {output}'

# AbSplice-DNA (no samples)
rule splicing_result_postprocess_absplice_dna_no_samples:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['absplice']['dna_no_samples'],
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['absplice']['dna_no_samples'],
    shell:
        'cp {input} {output}'

# AbSplice-RNA single cat
rule splicing_result_postprocess_absplice_rna_single_cat:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['absplice']['rna']['single_cat'],
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['absplice']['rna']['single_cat'],
    shell:
        'cp {input} {output}'

# AbSplice-RNA all cats
rule splicing_result_postprocess_absplice_rna_all_cats:
    input:
        model = OUTPUT_DIR_SPLICING_RAW + config_static['splicing_pred']['models']['absplice']['rna']['all_cats'],
    output:
        model_postprocess = OUTPUT_DIR_SPLICING_POSTPROCESS + config_static['splicing_pred']['models']['absplice']['rna']['all_cats'],
    shell:
        'cp {input} {output}'


list_outputs = list()
if 'mmsplice' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_mmsplice.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'mtsplice' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_mtsplice.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'spliceai' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_spliceai.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'cadd_splice' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_cadd_splice.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'mmsplice_splicemap' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_mmsplice_splicemap.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'mmsplice_splicemap_ref_psi' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_mmsplice_splicemap_ref_psi.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'spliceai_splicemap' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_spliceai_splicemap.output, 
        vcf_id=wildcard_vcf_id), # SpliceAI + SpliceMap
    )
if 'spliceai_splicemap_ref_psi' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_spliceai_splicemap_ref_psi.output, 
        vcf_id=wildcard_vcf_id), # SpliceAI + SpliceMap + PSI_ref
    )
if 'cat_infer_single' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_cat_infer_single_cat.output, 
        vcf_id=wildcard_vcf_id, tissue_cat=config['tissues_cat']), # CAT infer (single CAT)
    )
if 'cat_infer_all' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_cat_infer_all_cats.output, 
        vcf_id=wildcard_vcf_id), # CAT infer (all CAT)
    )
if 'absplice_dna' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_absplice_dna.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-DNA
    )
if 'absplice_dna_no_samples' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_absplice_dna_no_samples.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-DNA (input)
    )
if 'absplice_rna_single' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_absplice_rna_single_cat.output, 
        vcf_id=wildcard_vcf_id, tissue_cat=config['tissues_cat']), # AbSplice-RNA single CAT
    )
if 'absplice_rna_all' in config['models']:
    list_outputs.append(
        expand(rules.splicing_result_postprocess_absplice_rna_all_cats.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-RNA all CATs
    )


rule all_splicing_result_postprocess:
    input:
        list_outputs
        
        
del OUTPUT_DIR_SPLICING_RAW
del OUTPUT_DIR_SPLICING_POSTPROCESS
del var_samples_df