import pandas as pd
from tqdm import tqdm

df = pd.concat(
    [pd.read_csv(i) for i in tqdm(snakemake.input['rare_variants'])]
)
df = df.set_index(['junctions', 'variant'])

_df = df.groupby(df.index).size()

df = df.loc[_df[_df <= snakemake.params['max_num_samples']].index]
df.to_csv(snakemake.output['private_variants'])

