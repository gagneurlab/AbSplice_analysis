import pandas as pd
from absplice.utils import get_abs_max_rows
from absplice_scripts.utils.mapping_utils import subset_tissues

df_pred = pd.read_parquet(snakemake.input['model']).reset_index()

if 'tissue_subset' in snakemake.params.keys() and snakemake.params['tissue_subset'] == True:
    tissue_map = pd.read_csv(snakemake.input['tissue_map'])
    tissue_map = dict(zip(tissue_map['tissue'], tissue_map['tissue_main']))
    df_pred = subset_tissues(df_pred, tissue_map, chosen_tissue=snakemake.wildcards['tissue_pred'])

# Aggregate predictions
groupby_index = snakemake.params['groupby_index']
df_pred_agg = get_abs_max_rows(
    df_pred.set_index(groupby_index), 
    groupby=groupby_index, 
    max_col=snakemake.params['max_col']
).reset_index()

df_pred_agg.to_parquet(snakemake.output['model_agg'])
