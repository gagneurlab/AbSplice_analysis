import pandas as pd
from tqdm import tqdm
import pandas as pd
from absplice.utils import get_abs_max_rows
from absplice_scripts.data.gtex_v8 import tissue_map_mtsplice, tissue_map

df_pred = pd.read_csv(snakemake.input['pred_mtsplice'])\
            .rename(columns={'ID': 'variant'})
df_vcf_annotation = pd.read_csv(snakemake.input['var_samples_df'])

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
df_pred = df_pred.set_index('variant').join(
    df_vcf_annotation.set_index('variant'), 
    how='inner')\
    .reset_index()

# Aggregate predictions
df_pred_agg = get_abs_max_rows(
    df_pred.set_index(['variant', 'gene_id', 'sample', 'tissue']), 
    groupby=['variant', 'gene_id', 'sample', 'tissue'], 
    max_col='delta_logit_psi_mtsplice'
)

df_pred_agg.to_csv(snakemake.output['result_mtsplice_variant'])
