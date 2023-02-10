import pandas as pd

df_universe = pd.read_csv(snakemake.input['universe'])
df_samples_with_cat = pd.read_csv(snakemake.input['samples_with_cat'])

df_universe = df_universe.set_index('sample').join(
    df_samples_with_cat.set_index('sample')).reset_index()
df_universe.to_csv(snakemake.output['universe_with_cat'], index=False)