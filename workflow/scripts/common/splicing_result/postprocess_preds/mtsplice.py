import pandas as pd
from absplice_scripts.data.gtex_v8 import tissue_map_mtsplice, tissue_map

df_pred = pd.read_parquet(snakemake.input['model'])\
            .rename(columns={'ID': 'variant'})

# clean up gene_ids
df_pred['gene_id'] = df_pred['gene_id'].apply(lambda x: x.split(';'))
df_pred = df_pred.explode('gene_id')
df_pred['gene_id'] = df_pred['gene_id'].apply(lambda x: x.split('.')[0])

# rename tissues
df_pred = df_pred.rename(columns=tissue_map_mtsplice)

# melt dataframe
index_cols = ['variant', 'region', 'gene_id', 'gene_name']
value_cols = [x for x in df_pred.columns if x in sorted(tissue_map.keys())]
df_pred = df_pred.melt(
    id_vars=index_cols, 
    value_vars=value_cols, 
    var_name='tissue', 
    value_name='delta_logit_psi_mtsplice'
)

# join sample info
if 'var_samples_df' in snakemake.input.keys():
    df_vcf_annotation = pd.read_csv(snakemake.input['var_samples_df'])

    df_pred = df_pred.set_index('variant').join(
        df_vcf_annotation.set_index('variant'), 
        how='inner')\
        .reset_index()

df_pred.to_parquet(
    snakemake.output['model_postprocess'], 
    index=False)
