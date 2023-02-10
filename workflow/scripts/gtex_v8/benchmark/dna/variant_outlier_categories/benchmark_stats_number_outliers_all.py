import pandas as pd

df = pd.concat(
    [pd.read_csv(i) for i in snakemake.input['benchmark_stats']]
)

replace_dict_outlier = snakemake.params['replace_dict_outliers']
replace_dict_variants = snakemake.params['replace_dict_variants']

df = df.replace(replace_dict_variants).replace(replace_dict_outlier)
df = df.sort_values(by='num_outliers', ascending=False)

df.to_csv(snakemake.output['benchmark_stats_all'], index=False)
