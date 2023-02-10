import pandas as pd
from absplice_scripts.data.prokisch_WGS import tissue_map

rule mito_generate_gene_tpm:
    input:
        gene_tpm_gtex = config_precomputed['gene_tpm_gtex'].format(gtex_version=config['gtex_version']),
    output:
        gene_tpm_wide = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm_wide'],
        gene_tpm = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm'],
    run:
        df = pd.read_csv(input.gene_tpm_gtex)
        df_gene_tpm = df[df['tissue'] == 'Cells_Cultured_fibroblasts']
        df_gene_tpm = df_gene_tpm.replace({'Cells_Cultured_fibroblasts': 'Prokisch_Fibroblasts'})
        df_gene_tpm.to_csv(output.gene_tpm, index=False)
        df_gene_tpm_wide = df_gene_tpm.drop(columns='tissue')
        df_gene_tpm_wide = df_gene_tpm_wide.rename(columns={'gene_tpm': 'Prokisch_Fibroblasts'})
        df_gene_tpm_wide.to_csv(output.gene_tpm_wide, index=False)
        
        
rule mito_generate_tissue_map:
    output:
        tissue_map = config['DROP']['tissue_map'],
    run:
        df = pd.DataFrame(tissue_map, index=[0]).transpose().reset_index().rename(columns={'index': 'tissue', 0: 'tissue_DROP'})
        df.to_csv(output.tissue_map, index=False)
        