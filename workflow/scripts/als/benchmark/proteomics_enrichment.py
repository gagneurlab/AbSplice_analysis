import pandas as pd
from tqdm import tqdm
import numpy as np

# # # Matching cutoffs
# df_selected_cutoffs = pd.read_csv(snakemake.input['absplice_matching_cutoffs'])

model_cutoffs = snakemake.params['model_cutoffs']

# # Read in data
df_proteomics = pd.read_csv(snakemake.input['proteomics'])

df_gtex = pd.concat(
    [pd.read_parquet(i) for i in tqdm(snakemake.input['benchmark'])]
)

# ## Join predictions with proteomics
index = ['gene_id', 'sample']

df_gtex = df_gtex.set_index(index).join(
    df_proteomics.set_index(index)).reset_index()

df_gtex = df_gtex[~df_gtex['zScore'].isna()]
df_gtex['outlier'] = df_gtex['outlier'].fillna(False)

df = df_gtex[[
    'gene_id', 'sample', 'outlier', 'log2fc', 'zScore', 'padjust',
    *[x for x in df_gtex.columns if x in model_cutoffs.keys()]
]]
df = df.rename(columns={'outlier': 'splice_outlier'})

# # Get enrichment
threshold = 0.05
df_plot = list()

# for model in model_cutoffs.keys():
for model in [x for x in df_gtex.columns if x in model_cutoffs.keys()]:
    for cutoff_cat, cutoff in model_cutoffs[model].items():
        _df = df.copy()
        _df['cutoff_category'] = cutoff_cat
        _df['cutoff'] = cutoff
        _df['above_cutoff'] = np.abs(_df[model]) > cutoff
        _df['num_preds_above_cutoff'] = _df['above_cutoff'].sum()
        _df['num_preds_below_cutoff'] = _df[_df['above_cutoff'] == False].shape[0]
        _df['num_splicing_outlier_above_cutoff'] = _df[_df['above_cutoff'] == True]['splice_outlier'].sum()
        _df['model'] = model
        df_plot.append(_df)
        
df_plot = pd.concat(df_plot)

df_plot['cutoff_category'] = df_plot['cutoff_category'].astype("category")
df_plot['cutoff_category'] = df_plot['cutoff_category'].cat.reorder_categories(['low', 'medium', 'high'])

df_plot.to_csv(snakemake.output['enrichment'], index=False)

# # Get proportion validated
zscore_proteomics_cutoff = -2
df_prot_true = pd.DataFrame(
    df_plot[
        df_plot['zScore'] <= zscore_proteomics_cutoff
    ].groupby(
        ['model', 'cutoff_category', 'above_cutoff']
    ).apply(
        lambda df: df.shape[0])
).reset_index().rename(columns={0: 'num_protein_low'})

df_prot_true = df_prot_true[
    df_prot_true['above_cutoff'] == True
]

df_stats = pd.DataFrame(
    df_plot.groupby(
        ['model', 'cutoff_category', 'above_cutoff']
    ).apply(
        lambda df: df.shape[0])
).reset_index().rename(columns={0: 'num'})

df_stats = df_stats[
    df_stats['above_cutoff'] == True
]

df_stats = df_stats.set_index(['model', 'cutoff_category', 'above_cutoff']).join(
    df_prot_true.set_index(['model', 'cutoff_category', 'above_cutoff'])).reset_index()

df_stats['proportion_true'] = df_stats['num_protein_low']/df_stats['num']

df_all_genes = pd.DataFrame({
    'model': ['All genes'],
    'cutoff_category': ['high'],
    'above_cutoff': [True],
    'num': [df_plot.shape[0]],
    'num_protein_low': [df_plot[df_plot['zScore'] <= -2].shape[0]],
})
df_all_genes['proportion_true'] = df_all_genes['num_protein_low']/df_all_genes['num']
df_all_genes

df_stats = pd.concat([df_stats, df_all_genes])

df_stats.to_csv(snakemake.output['true_preds'], index=False)
