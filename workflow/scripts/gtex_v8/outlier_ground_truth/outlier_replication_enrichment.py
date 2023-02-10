import pandas as pd

df_var_outlier_dist = pd.read_csv(snakemake.input['variant_outlier_dist'])
df_var_outlier_dist = df_var_outlier_dist[['junctions', 'sample', 'variant', 'gene_id', 'abs_Distance']]

df_var_outlier_dist = df_var_outlier_dist[
    (df_var_outlier_dist['abs_Distance'] <= int(snakemake.wildcards['rare_var_dist']))
]

df_outlier_multiple_tissues = pd.read_csv(snakemake.input['outlier_multiple_tissues_junction_level'])

df_outlier_multiple_tissues[
    df_outlier_multiple_tissues['tissue'] == snakemake.wildcards['tissue']
]

index = ['junctions', 'sample', 'gene_id']

df_all = df_var_outlier_dist.set_index(index).join(
    df_outlier_multiple_tissues.set_index(index)[['num_tissues_per_outlier']])

df_all['num_tissues_per_outlier'] = df_all['num_tissues_per_outlier'].fillna(0)
df_all['num_tissues_per_outlier'] = df_all['num_tissues_per_outlier'].astype(int)

df_all = df_all.reset_index().drop_duplicates()
df_all.to_csv(snakemake.output['variant_outlier_dist_multiple_tissues'], index=False)
