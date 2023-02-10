import numpy as np
import pandas as pd
from tqdm import tqdm
from splicemap import SpliceCountTable as CountTable

ct = CountTable.read_csv(snakemake.input['count_table_updated'], 
                         name=snakemake.params['tissue'])
df_private_variants = pd.read_csv(snakemake.input['private_variants'])

junctions = set(ct.junctions)
samples = set(ct.samples)

df = pd.DataFrame(
    np.zeros((len(junctions), len(samples))),
    columns=samples,
    index=junctions
)

for _, row in tqdm(df_private_variants.iterrows()):
    if row['junctions'] in junctions and row['sample'] in samples:
        df.loc[row['junctions'], row['sample']] = 1

df.to_csv(snakemake.output['train_test_split'])
