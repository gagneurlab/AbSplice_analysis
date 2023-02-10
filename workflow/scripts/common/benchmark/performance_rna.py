import pandas as pd
from tqdm import tqdm
import pickle
import os
import re
from absplice_scripts.visualization.benchmark import get_performance
from absplice_scripts.utils.mapping_utils import subset_tissues

valid_cols = [
    'outlier',
    *snakemake.params['model_dict'].values(),
]

def subset_samples_paired_with_cat(_df):
    tissues_cat = snakemake.wildcards['cat_pairing'].split('__')
    for tissue_cat in tqdm(tissues_cat):
        _df = _df[_df[tissue_cat] == True]
    return _df

def median_cutoff(_df):
    # median_n cutoff for SpliceMap results
    _df.loc[
        _df['median_n_mmsplice_splicemap'] < snakemake.params['median_n_cutoff'], 
        'delta_logit_psi_mmsplice_splicemap'] = 0
    _df.loc[
        _df['median_n_mmsplice_splicemap_ref_psi'] < snakemake.params['median_n_cutoff'], 
        'delta_psi_mmsplice_splicemap_ref_psi'] = 0
    
    _df.loc[
        (_df['psi3_median_n_spliceai_splicemap'] < snakemake.params['median_n_cutoff']) \
        | (_df['psi5_median_n_spliceai_splicemap'] < snakemake.params['median_n_cutoff']),
        'delta_score_spliceai_splicemap'] = 0
    _df.loc[
        (_df['psi3_median_n_spliceai_splicemap_ref_psi'] < snakemake.params['median_n_cutoff']) \
        | (_df['psi5_median_n_spliceai_splicemap_ref_psi'] < snakemake.params['median_n_cutoff']),
        'delta_score_spliceai_splicemap_ref_psi'] = 0
    
    _df.loc[
        _df['median_n_cat_mmsplice_splicemap_cat'] < snakemake.params['median_n_cutoff'], 
        'delta_psi_cat_mmsplice_splicemap_cat'] = 0
    return _df

df_benchmark = list()
for i in tqdm(snakemake.input['benchmark']):
    _df = pd.read_parquet(i)
    _df = subset_samples_paired_with_cat(_df)
    _df = median_cutoff(_df)
    df_benchmark.append(_df[valid_cols])
df_benchmark = pd.concat(df_benchmark)

df_performance, performance = get_performance(
    df = df_benchmark, 
    outlier_column = snakemake.params['outlier_column'], 
    model_dict = snakemake.params['model_dict'],
)

df_performance.to_csv(snakemake.output['df_performance'], index=False)
pickle.dump(performance, open(snakemake.output['aps_performance'], 'wb'))
