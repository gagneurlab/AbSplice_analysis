# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python [conda env:anaconda-absplice_paper]
#     language: python
#     name: conda-env-anaconda-absplice_paper-py
# ---

# %%
from absplice_scripts.dataclasses.junction import get_splice_site_intervals, \
    get_unique_splice_site_intervals_in_event, intervals_to_pyranges
from splicemap import SpliceCountTable as CountTable
import pandas as pd
from typing import List, Union, Iterable, Iterator
import pyranges as pr
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors import variants_to_pyranges
from tqdm import tqdm

from splicemap.count_table import SpliceCountTable

# %%
try:
    snakemake
except NameError:
    import os
#     os.chdir('../..')
    os.chdir('/data/nasif12/home_if12/wagnern/Projects/gitlab_gagneurlab/absplice/workflow/gtex_v8/from_general_workflow/outlier_ground_truth_all')

    from snakemk_util import load_rule_args

    snakemake = load_rule_args(
        snakefile = 'Snakefile',
        rule_name = 'outlier_var_nearest_variant_centric',
        root=os.getcwd(),
        default_wildcards={
#             'features': ['delta_score', 'delta_psi'],
            # 'tissue': 'Heart_Left_Ventricle',
            # 'vcf_id': 'chr14',
            # 'delta_psi_cutoff': 0.1,
            # 'outlier_type': 'psi5__psi3',
            
            'tissue': 'Adipose_Subcutaneous',
            'vcf_id': 'chr1',
            'delta_psi_cutoff': 0.3,
            'outlier_type': 'psi5__psi3__theta',
        }
    )

# %%
snakemake.output

# %% [markdown]
# # Outliers

# %%
snakemake.input['outliers_signif']

# %%
df_outlier = pd.read_csv(snakemake.input['outliers_signif'])

df_outlier['junctions'] = df_outlier['junctions_j']
df_outlier['events'] = df_outlier['events_j']

df_outlier = df_outlier.set_index('junctions')
df_outlier['events'] = df_outlier['events'].str.split(';')

# %%
df_outlier = df_outlier[
    df_outlier['seqnames_j'] == snakemake.wildcards['vcf_id']
]

# %%
sorted(set(df_outlier['seqnames_j']))

# %%
df_outlier.shape

# %% [markdown]
# # Rare variants

# %%
df_rare_vars = pd.read_csv(snakemake.input['rare_vars'])

df_rare_vars.shape

# %%
common_samples = sorted(set(df_outlier['sample']).intersection(
    set(df_rare_vars['sample'])))

# %%
len(common_samples)


# %% [markdown]
# # Distance junctions to nearest variant

# %%
def junc_var_dist(df_outlier, df_rare_vars, sample):
    df_outlier = df_outlier[df_outlier['sample'] == sample]
    df_rare_vars = df_rare_vars[df_rare_vars['sample'] == sample]
    
    # get unique variants into pyranges
    v_list = df_rare_vars['variant'].apply(lambda x: Variant.from_str(x)).values
    v_list = list(set(v_list))
    pr_rare_vars = variants_to_pyranges(v_list)
    
    # get unique intervals of outlier events
    df_outlier['interval'] = df_outlier['events'].apply(lambda x: get_unique_splice_site_intervals_in_event(x, overhang=(0, 0)))
    df_intervals_all = df_outlier.explode('interval')[['interval']]
        
    # get unique intervals into pyranges
    i_list = df_intervals_all['interval'].values
    i_list = list(set(i_list))
    pr_intervals_unique = intervals_to_pyranges(i_list)
    df_intervals_unique = pr_intervals_unique.df
    
    # get nearest distance of junctions to rare variants
    # NOTE: HERE VARIANT COMES FIRST, ALL VARIANTS HAVE TO BE ANNOTATED
    pr_intervals_with_rare = pr_rare_vars.k_nearest(pr_intervals_unique, k=1)
    df_intervals_with_rare = pr_intervals_with_rare.df
    
    # interval to string for joining
    df_intervals_with_rare['variant'] = df_intervals_with_rare['variant'].apply(lambda x: x.__str__())
    df_intervals_all['interval'] = df_intervals_all['interval'].apply(lambda x: x.__str__())
    df_intervals_unique['interval'] = df_intervals_unique['interval'].apply(lambda x: x.__str__())
    df_intervals_with_rare['interval'] = df_intervals_with_rare['interval'].apply(lambda x: x.__str__())
    
    # get junction information from df_intervals_all
    df_intervals_all = df_intervals_all.reset_index().set_index('interval')
    df_intervals_with_rare = df_intervals_with_rare.set_index('interval')
    df_junctions_with_rare_event = df_intervals_all.join(
        df_intervals_with_rare, how='inner').reset_index()
    
    # get variant informations from df_rare_vars
    df_junctions_with_rare_event['variant'] = df_junctions_with_rare_event['variant'].apply(lambda x: x.__str__())
    df_junctions_with_rare_event = df_junctions_with_rare_event.set_index('variant').join(
        df_rare_vars.set_index('variant')).reset_index()
        
    return df_junctions_with_rare_event


# %%
df_junc_vars = list()
for sample in tqdm(common_samples):
    df_junc_vars.append(junc_var_dist(df_outlier, df_rare_vars, sample))
    assert len(set(df_rare_vars[df_rare_vars['sample'] == sample].reset_index()['variant']).difference(
        set(df_junc_vars[-1]['variant']))) == 0
df_junc_vars = pd.concat(df_junc_vars)


# %% [markdown]
# # get minimum distance of variant and junction

# %%
def get_abs_min_rows(df, groupby, min_col, dropna=True):
    df = df.reset_index()
    _df = df.copy()
    _df[min_col] = _df[min_col].abs()
    min_scores = _df.groupby(groupby, dropna=dropna)[min_col].idxmin()
    return df.iloc[min_scores.values].set_index(groupby)


# %%
# get minimum distance of variant and junction
df_junc_vars_min = get_abs_min_rows(
    df_junc_vars.set_index(['sample', 'variant']),
    min_col='Distance',
    groupby=['sample', 'variant']
).reset_index()

# %%
df_junc_vars_min['abs_Distance'] = abs(df_junc_vars_min['Distance'])

# %% [markdown]
# # annotate gene

# %%
df_junc_vars_min = df_junc_vars_min.set_index(['junctions', 'sample']).join(
    df_outlier.reset_index()[['junctions', 'sample', 'gene_id', 'pValueGene_g']].drop_duplicates().set_index(['junctions', 'sample'])).reset_index()

# %% [markdown] tags=[]
# # Save results

# %%
snakemake.output

# %%
df_junc_vars_min.to_csv(snakemake.output['var_junc_nearest'], index=False)
