import pandas as pd

df_pred = pd.read_parquet(snakemake.input['model'])

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
