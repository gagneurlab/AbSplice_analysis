OUTPUT_DIR_RAW = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
OUTPUT_DIR_GENE = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['gene_level']

import pandas as pd            
# --------------------------------gene level-------------------------------
# AbSplice
rule splicing_result_gene_level_absplice_dna_gtex_splicemaps:
    input:
        pred_absplice_dna = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_absplice_dna_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['dna'],
    script:
        "../gene_level/gene_absplice_dna.py"
               
            
rule splicing_result_gene_level_absplice_rna_gtex_splicemaps:
    input:
        pred_absplice_rna = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['single_cat'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_absplice_rna_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['single_cat'],
    script:
        "../gene_level/gene_absplice_rna.py"
        

rule splicing_result_gene_level_absplice_rna_gtex_splicemaps_all_cats:
    input:
        pred_absplice_rna = OUTPUT_DIR_RAW + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['all_cats'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_absplice_rna_gene = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['absplice']['gtex_splicemaps']['rna']['all_cats'],
    script:
        "../gene_level/gene_absplice_rna.py"