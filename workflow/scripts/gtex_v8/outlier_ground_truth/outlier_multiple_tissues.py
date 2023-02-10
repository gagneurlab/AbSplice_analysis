from tqdm import tqdm
import pandas as pd

# sample annotation, number tissues per sample
df_tissue_map = pd.read_csv(snakemake.input['tissue_map'])
tissue_map = dict(zip(df_tissue_map['tissue_DROP'], df_tissue_map['tissue']))

df_anno = pd.read_csv(snakemake.input['sample_anno'], sep='\t')
df_anno['DROP_GROUP'] = df_anno['DROP_GROUP'].map(tissue_map)

df_anno = df_anno[
    (~df_anno['INDIVIDUAL_ID'].isna())
    & (~df_anno['RNA_ID'].isna())
]
df_anno = df_anno[['INDIVIDUAL_ID', 'DROP_GROUP']]\
            .rename(columns={'DROP_GROUP': 'tissue', 'INDIVIDUAL_ID': 'sample'})
df_anno = df_anno.dropna()

df_tissues_sample = df_anno.groupby('sample').apply(lambda x: sorted(list(set(x['tissue']))))
df_tissues_sample = pd.DataFrame(df_tissues_sample).rename(columns={0: 'tissues_per_sample'})
df_tissues_sample['num_tissues_per_sample'] = df_tissues_sample['tissues_per_sample'].apply(lambda x: len(x))

valid_samples = sorted(list(set(df_tissues_sample[df_tissues_sample['num_tissues_per_sample'] > 1].index)))

# Outliers
df_outliers_signif = pd.concat(
    [pd.read_csv(i).assign(tissue=t) for i, t in tqdm(zip(snakemake.input['outliers_signif'], snakemake.params['tissues']))]
)

df_outliers_signif_junc = df_outliers_signif[['gene_id', 'junctions_j', 'sample', 'tissue']].drop_duplicates()
df_outliers_signif_junc = df_outliers_signif_junc.rename(columns={'junctions_j': 'junctions'})

# Gene level
df_outliers_signif = df_outliers_signif[['gene_id', 'sample', 'tissue']].drop_duplicates()

df_all_tissue_outlier = df_outliers_signif.groupby(['sample', 'gene_id']).apply(lambda x: sorted(list(set(x['tissue']))))
df_all_tissue_outlier = pd.DataFrame(df_all_tissue_outlier).rename(columns={0: 'tissues_per_outlier'})
df_all_tissue_outlier['num_tissues_per_outlier'] = df_all_tissue_outlier['tissues_per_outlier'].apply(lambda x: len(x))
df_all_tissue_outlier = df_all_tissue_outlier[df_all_tissue_outlier.index.get_level_values('sample').isin(valid_samples)]

df_multiple_tissue_outlier = df_all_tissue_outlier[df_all_tissue_outlier['num_tissues_per_outlier'] > 1]
df_multiple_tissue_outlier['tissue'] = df_multiple_tissue_outlier['tissues_per_outlier']
df_multiple_tissue_outlier = df_multiple_tissue_outlier.explode('tissue')
df_multiple_tissue_outlier.to_csv(snakemake.output['outlier_multiple_tissues'])

# Junction level
df_all_tissue_outlier_junc = df_outliers_signif_junc.groupby(['sample', 'gene_id', 'junctions']).apply(lambda x: sorted(list(set(x['tissue']))))
df_all_tissue_outlier_junc = pd.DataFrame(df_all_tissue_outlier_junc).rename(columns={0: 'tissues_per_outlier'})
df_all_tissue_outlier_junc['num_tissues_per_outlier'] = df_all_tissue_outlier_junc['tissues_per_outlier'].apply(lambda x: len(x))
df_all_tissue_outlier_junc = df_all_tissue_outlier_junc[df_all_tissue_outlier_junc.index.get_level_values('sample').isin(valid_samples)]

df_multiple_tissue_outlier_junc = df_all_tissue_outlier_junc[df_all_tissue_outlier_junc['num_tissues_per_outlier'] > 1]
df_multiple_tissue_outlier_junc['tissue'] = df_multiple_tissue_outlier_junc['tissues_per_outlier']
df_multiple_tissue_outlier_junc = df_multiple_tissue_outlier_junc.explode('tissue')
df_multiple_tissue_outlier_junc.to_csv(snakemake.output['outlier_multiple_tissues_junction_level'])