import pandas as pd
from tqdm import tqdm

df_tissue_map = pd.read_csv(snakemake.input['tissue_map'])
tissue_map = dict(zip(df_tissue_map['tissue_DROP'], df_tissue_map['tissue']))

df_anno = pd.read_csv(snakemake.input['sample_anno'], sep='\t')
df_anno['DROP_GROUP'] = df_anno['DROP_GROUP'].map(tissue_map)
df_anno = df_anno[['INDIVIDUAL_ID', 'DROP_GROUP']]\
    .rename(columns={'INDIVIDUAL_ID': 'sample', 'DROP_GROUP': 'tissue'})

df_tissues_per_sample = pd.DataFrame(df_anno.groupby('sample').apply(
    lambda df: sorted(set(df['tissue']))))\
    .rename(columns={0:'tissues_per_sample'})\
    .reset_index()

for tissue_cat in tqdm(snakemake.params['tissues_cat']):
    df_tissues_per_sample.loc[:, tissue_cat] = df_tissues_per_sample.apply(
        lambda df: tissue_cat in df['tissues_per_sample'], axis=1)
    
df_tissues_per_sample = df_tissues_per_sample.drop(columns='tissues_per_sample')
df_tissues_per_sample.to_csv(snakemake.output['samples_with_cat'], index=False)