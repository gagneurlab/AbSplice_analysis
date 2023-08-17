OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
ABSPLICE_INPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['combined_input']

##==================================SPLICING PREDICTIONS (BASELINE MODELS)=========================================    
# MMSplice
rule splicing_pred_mmsplice:
    input:
        vcf = config['vcf'],
        fasta = config['fasta'],
        gtf = config['gtf'],
    params:
        genome = config['genome'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1
    output:
        result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice'])
    script:
        "../mmsplice_exon.py"
        
# MTSplice
rule splicing_pred_mtsplice:
    input:
        vcf = config['vcf'],
        fasta = config['fasta'],
        gtf = config['gtf'],
#     conda: 'env/mtsplice_absplice.yaml' #mtsplice_absplice
    params:
        genome = config['genome'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1
    output:
        result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['mtsplice'])
    script:
        "../mtsplice.py"

# SpliceAI
rule splicing_pred_spliceai:
    input:
        vcf = config['vcf'],
        db = config_precomputed['spliceai']['db'].format(
            genome=config['genome']),
        fasta = config['fasta'],
    # conda:
    #     '/opt/modules/i12g/anaconda/envs/spliceai_gpu'
    params:
        lookup_only = config['spliceai']['lookup_only'],
        genome = config['genome'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1,
        gpu = 1,
    output:
        result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'])
    script:
        "../spliceai.py"
        
# CADD-Splice        
genome_mapper = {
    'hg38': 'GRCh38',
    'hg19': 'GRCh37',
}

rule splicing_pred_cadd_splice:
    input:
        vcf = config['vcf'],
    conda: '/opt/modules/i12g/anaconda/envs/cadd_splice'
    params:
        genome = genome_mapper[config['genome']],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1
    output:
        result = OUTPUT_DIR + config_static['splicing_pred']['models']['cadd_splice']
    shell:
        '''/data/nasif12/home_if12/wagnern/Projects/CADD-scripts/CADD.sh -a -g "{params.genome}" -o {output.result} {input.vcf}'''
    

##==================================SPLICING PREDICTIONS (WITH SPLICEMAP)=========================================    
# MMSplice + SpliceMap
rule splicing_pred_mmsplice_splicemap:
    input:
        vcf = config['vcf'],
        fasta = config['fasta'],
        splicemap_5 = expand(config_static['junction_annotation']['splicemap']['psi5'],
                             genome=config['genome'], tissue=config['tissues']),
        splicemap_3 = expand(config_static['junction_annotation']['splicemap']['psi3'],
                             genome=config['genome'], tissue=config['tissues']),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 4
    output:
        result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap'])
    script:
        "../mmsplice_splicemap.py"
        
# SpliceAI + SpliceMap (single tissue)
rule splicing_pred_spliceai_splicemap_tissue:
    input:
        spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
        splicemap_5 = expand(config_static['junction_annotation']['splicemap']['psi5'],
                             genome=config['genome'], tissue='{tissue}'),
        splicemap_3 = expand(config_static['junction_annotation']['splicemap']['psi3'],
                             genome=config['genome'], tissue='{tissue}'),
        gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
    params:
        tissue = '{tissue}'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 5000,
        threads = 1,
    output:
        spliceai_splicemap = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['tissue'])
    script:
        "../spliceai_splicemap.py"

# SpliceAI + SpliceMap (all tissues)      
rule splicing_pred_spliceai_splicemap_all:
    input:
        df_tissue = expand(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['tissue'],
                           vcf_id='{vcf_id}', tissue=config['tissues']),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1,
    output:
        df_all = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['all'])
    run:
        import pandas as pd
        from tqdm import tqdm
        df = pd.concat([pd.read_parquet(i) for i in tqdm(input.df_tissue)])
        df.to_parquet(output.df_all, index=False)
        
# SpliceAI + SpliceMap + PSI_ref (single tissue)    
rule splicing_pred_spliceai_splicemap_ref_psi_tissue:
    input:
        spliceai_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['tissue']
    params:
        tissue = '{tissue}',
        ref_psi_tolerance = 0.05
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1,
    output:
        spliceai_splicemap_ref_psi = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['tissue'])
    script:
        "../spliceai_splicemap_ref_psi.py"
        
# SpliceAI + SpliceMap + PSI_ref (all tissues)         
rule splicing_pred_spliceai_splicemap_ref_psi_all:
    input:
        df_tissue = expand(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['tissue'],
                           vcf_id='{vcf_id}', tissue=config['tissues']),
    output:
        df_all = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['all'])
    run:
        import pandas as pd
        from tqdm import tqdm
        df = pd.concat([pd.read_parquet(i) for i in tqdm(input.df_tissue)])
        df.to_parquet(output.df_all, index=False)

# CAT infer (single CAT)     
if 'cat_infer_single' in config['models']:       
    rule splicing_pred_infer_cat:
        input:
            var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
            mmsplice_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap'],
            splicemap_5 = expand(config_static['junction_annotation']['splicemap']['psi5'],
                                genome=config['genome'], tissue=config['tissues']),
            splicemap_3 = expand(config_static['junction_annotation']['splicemap']['psi3'],
                                genome=config['genome'], tissue=config['tissues']),
            cat_count_table = expand(OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['updated'],
                                    tissue='{tissue_cat}', genome=config['genome'])[0],
        params:
            tissue_cat = '{tissue_cat}',
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 64000,
            threads = 1
        output:
            result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat'])
        script:
            "../infer_cat.py"
        
# CAT infer (all CATs)
if 'cat_infer_all' in config['models']:
    rule splicing_pred_mmsplice_splicemap_cat_combine:
        input:
            pred_mmsplice_cat = expand(OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat'],
                                    tissue_cat=config['tissues_cat'], vcf_id='{vcf_id}'),
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 8000,
            threads = 1
        output:
            result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['all_cats'])
        run:
            df_mmsplice_cat = pd.concat([pd.read_parquet(i).reset_index() for i in input.pred_mmsplice_cat])
            df_mmsplice_cat.to_parquet(output.result, partition_cols='tissue_cat', index=False)
        
        
# --------------------------------AbSplice input------------------------------- 
# AbSplice-DNA (input)
if 'absplice_dna' in config['models'] or 'absplice_dna_input' in config['models']:
    rule splicing_pred_absplice_dna_input:
        input:
            pred_mmsplice = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap'],
            pred_spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
            gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
            gene_tpm = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm'],
            var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 64000,
            threads = 1
        output:
            absplice_dna_input = directory(ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna'])
        script:
            "../absplice_dna_input.py"
        
# AbSplice-RNA single CAT (input)    
if 'absplice_rna_single' in config['models'] or 'absplice_rna_single_input' in config['models']:
    rule splicing_pred_absplice_rna_single_cat_input:
        input:
            absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna'],
            pred_mmsplice_cat = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat'],
            CAT_pval = expand(OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['tissue_cat_pval'],
                            tissue='{tissue_cat}', vcf_id='{vcf_id}'),
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 64000,
            threads = 1
        output:
            absplice_rna_input = directory(ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['single_cat'])
        script:
            "../absplice_rna_input.py"
         
# AbSplice-RNA all CATs (input)          
if 'absplice_rna_all' in config['models'] or 'absplice_rna_all_input' in config['models']:
    rule splicing_pred_absplice_rna_all_cats_input:
        input:
            absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna'],
            pred_mmsplice_cat = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['all_cats'],
            CAT_pval = expand(OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['tissue_cat_pval'],
                            tissue=config['tissues_cat'], vcf_id='{vcf_id}'),
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 64000,
            threads = 1
        output:
            absplice_rna_input = directory(ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['all_cats'])
        script:
            "../absplice_rna_input.py"


# --------------------------------predict AbSplice------------------------------- 
if 'absplice_dna' in config['models']:
    rule splicing_pred_absplice_dna:
        input:
            absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna'],
            absplice_model = config_precomputed['absplice']['dna'],
        params:
            median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
            tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 32000,
            threads = 1
        output:
            absplice_dna_pred = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna'])
        script:
            "../absplice_dna.py"
        
        
if 'absplice_rna_single' in config['models']:
    rule splicing_pred_absplice_rna_single_cat:
        input:
            absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['single_cat'],
            absplice_model = config_precomputed['absplice']['rna'],
        params:
            median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
            tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 32000,
            threads = 1
        output:
            absplice_rna_pred = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['single_cat'])
        script:
            "../absplice_rna.py"
        

if 'absplice_rna_all' in config['models']:        
    rule splicing_pred_absplice_rna_all_cats:
        input:
            absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['all_cats'],
            absplice_model = config_precomputed['absplice']['rna'],
        params:
            median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
            tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 32000,
            threads = 1
        output:
            absplice_rna_pred = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['all_cats'])
        script:
            "../absplice_rna.py"


# --------------------------------AbSplice-DNA (no samples))-------------------------------     
if config['use_gene_tpm'] == True:
    gene_tpm_used = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm'],
else:
    gene_tpm_used = []

if 'absplice_dna_no_samples' in config['models']: 
    rule splicing_pred_absplice_dna_input_no_samples:
        input:
            pred_mmsplice = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap'],
            pred_spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
            gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
            gene_tpm = gene_tpm_used,
            splicemap_5 = expand(config_static['junction_annotation']['splicemap']['psi5'],
                                genome=config['genome'], tissue=config['tissues']),
            splicemap_3 = expand(config_static['junction_annotation']['splicemap']['psi3'],
                                genome=config['genome'], tissue=config['tissues']),
        params:
            tissues = config['tissues']
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 48000,
            threads = 1
        output:
            absplice_dna_input = directory(ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna_no_samples'])
        script:
            "../absplice_dna_input_no_samples.py"
        

if 'absplice_dna_no_samples' in config['models']: 
    rule splicing_pred_absplice_dna_no_samples:
        input:
            absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna_no_samples'],
            absplice_model = config_precomputed['absplice']['dna'],
        params:
            median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
            tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 48000,
            threads = 1
        output:
            absplice_dna_pred = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna_no_samples'])
        script:
            "../absplice_dna.py"


list_outputs = list()
if 'mmsplice' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_mmsplice.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'mtsplice' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_mtsplice.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'spliceai' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_spliceai.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'cadd_splice' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_cadd_splice.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'mmsplice_splicemap' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_mmsplice_splicemap.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap
    )
if 'mmsplice_splicemap_ref_psi' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_mmsplice_splicemap.output, 
        vcf_id=wildcard_vcf_id), # MMSplice + SpliceMap + PSI_ref
    )
if 'spliceai_splicemap' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_spliceai_splicemap_all.output, 
        vcf_id=wildcard_vcf_id), # SpliceAI + SpliceMap
    )
if 'spliceai_splicemap_ref_psi' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_spliceai_splicemap_ref_psi_all.output, 
        vcf_id=wildcard_vcf_id), # SpliceAI + SpliceMap + PSI_ref
    )
if 'cat_infer_single' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_infer_cat.output, 
        vcf_id=wildcard_vcf_id, tissue_cat=config['tissues_cat']), # CAT infer (single CAT)
    )
if 'cat_infer_all' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_mmsplice_splicemap_cat_combine.output, 
        vcf_id=wildcard_vcf_id), # CAT infer (all CAT)
    )
if 'absplice_dna__input' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_dna_input.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-DNA (input)
    )
if 'absplice_dna_no_samples__input' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_dna_input_no_samples.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-DNA (input)
    )
if 'absplice_rna_single__input' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_rna_single_cat_input.output, 
        vcf_id=wildcard_vcf_id, tissue_cat=config['tissues_cat']), # AbSplice-RNA single CAT (input)
    )
if 'absplice_rna_all__input' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_rna_all_cats_input.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-RNA all CATs (input)
    )
if 'absplice_dna' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_dna.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-DNA
    )
if 'absplice_dna_no_samples' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_dna_no_samples.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-DNA (input)
    )
if 'absplice_rna_single' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_rna_single_cat.output, 
        vcf_id=wildcard_vcf_id, tissue_cat=config['tissues_cat']), # AbSplice-RNA single CAT
    )
if 'absplice_rna_all' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_rna_all_cats.output, 
        vcf_id=wildcard_vcf_id), # AbSplice-RNA all CATs
    )


rule all_splicing_pred:
    input:
        list_outputs
        
        
del OUTPUT_DIR
del ABSPLICE_INPUT_DIR