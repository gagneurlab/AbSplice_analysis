import pandas as pd
from tqdm import tqdm

# variant outlier junction distance
variant_outlier_dist = pd.read_csv(snakemake.input['variant_outlier_dist'])

variant_outlier_dist = variant_outlier_dist[
    variant_outlier_dist['abs_Distance'] <= 250]

variant_outlier_dist = variant_outlier_dist.set_index(['junctions', 'sample', 'gene_id'])

# outlier outcome
junction_all_annotated = pd.read_csv(snakemake.input['outlier_outcome']).set_index(['junctions', 'sample', 'gene_id'])

# join
outliers = variant_outlier_dist.join(
            junction_all_annotated[['annotatedJunction', 'blacklist', 'aberrantSpliceType', 'causesFrameshift']]).reset_index()


annotation_cols = ['annotatedJunction', 'blacklist', 'aberrantSpliceType', 'causesFrameshift', 'junctions']
index_cols = ['variant', 'gene_id', 'sample']

df_all_annotations = list()
for col in tqdm(annotation_cols):
    _df = pd.DataFrame(outliers.groupby(index_cols).apply(lambda df: list(set(df[col])))).rename(columns={0: col})
    df_all_annotations.append(_df)
df_all_annotations = pd.concat(df_all_annotations, axis=1)

df_all_annotations = df_all_annotations.fillna('')

if df_all_annotations.shape[0] == 0:
    df_all_annotations = pd.DataFrame(columns=[*index_cols, *[x for x in annotation_cols if x != 'junctions']]).set_index(index_cols)

def get_abs_min_rows(df, groupby, min_col, dropna=True):
    df = df.reset_index()
    _df = df.copy()
    _df[min_col] = _df[min_col].abs()
    min_scores = _df.groupby(groupby, dropna=dropna)[min_col].idxmin()
    return df.iloc[min_scores.values].set_index(groupby)

# get minimum distance of variant and junction
df_min_dist = get_abs_min_rows(
    outliers.set_index(['variant', 'gene_id', 'sample']),
    min_col='Distance',
    groupby=['variant', 'gene_id', 'sample']
)

df_min_dist = df_min_dist[['Distance', 'abs_Distance']]

df_all = df_all_annotations.join(df_min_dist)
df_all.reset_index().to_csv(snakemake.output['outliers_with_annotation'], index=False)