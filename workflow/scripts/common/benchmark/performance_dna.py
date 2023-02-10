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
    return _df

df_benchmark = list()
for i in tqdm(snakemake.input['benchmark']):
    _df = pd.read_parquet(i)
    _df = median_cutoff(_df)
    df_benchmark.append(_df[valid_cols])
df_benchmark = pd.concat(df_benchmark)

# # subset tissue for single tissue performance
# if snakemake.params['subset_tissue'] == True:
#     tissue_map = pd.read_csv(snakemake.input['tissue_map'])
#     tissue_map = dict(zip(tissue_map['tissue'], tissue_map['tissue_main']))
#     df_benchmark = subset_tissues(
#         df_benchmark, tissue_map, chosen_tissue=snakemake.wildcards['tissue'])
#     assert len(set(df_benchmark['tissue']).difference(set([snakemake.wildcards['tissue']]))) == 0

df_performance, performance = get_performance(
    df = df_benchmark, 
    outlier_column = snakemake.params['outlier_column'], 
    model_dict = snakemake.params['model_dict'],
)

df_performance.to_csv(snakemake.output['df_performance'], index=False)
pickle.dump(performance, open(snakemake.output['aps_performance'], 'wb'))