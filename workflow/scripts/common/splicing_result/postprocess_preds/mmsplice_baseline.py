import pandas as pd
from absplice.utils import get_abs_max_rows

df_pred = pd.read_parquet(snakemake.input['model'])\
            .rename(columns={'ID': 'variant'})
     
# clean up gene_ids
df_pred['gene_id'] = df_pred['gene_id'].apply(lambda x: x.split(';'))
df_pred = df_pred.explode('gene_id')
df_pred['gene_id'] = df_pred['gene_id'].apply(lambda x: x.split('.')[0])
       
# add sample info
if 'var_samples_df' in snakemake.input.keys():
    df_vcf_annotation = pd.read_csv(snakemake.input['var_samples_df'])
        
    df_pred = df_pred.set_index('variant').join(
        df_vcf_annotation.set_index('variant'), how='inner')\
        .reset_index()

df_pred.to_parquet(snakemake.output['model_postprocess'], index=False)
