import pandas as pd
from tqdm import tqdm
import numpy as np

def median_n_filter(_df):
    # median_n cutoff for SpliceMap results
    _df.loc[
        _df['median_n_mmsplice_splicemap'] < snakemake.params['median_n_cutoff'], 
        'delta_logit_psi_mmsplice_splicemap'] = 0
    _df.loc[
        _df['median_n_mmsplice_splicemap'] < snakemake.params['median_n_cutoff'], 
        'delta_psi_mmsplice_splicemap'] = 0

    _df.loc[
        (_df['psi3_median_n_spliceai_splicemap'] < snakemake.params['median_n_cutoff']) | (_df['psi5_median_n_spliceai_splicemap'] < snakemake.params['median_n_cutoff']),
        'delta_score_spliceai_splicemap'] = 0
    _df.loc[
        (_df['psi3_median_n_spliceai_splicemap_ref_psi'] < snakemake.params['median_n_cutoff']) | (_df['psi5_median_n_spliceai_splicemap_ref_psi'] < snakemake.params['median_n_cutoff']),
        'delta_score_spliceai_splicemap_ref_psi'] = 0
    
    return _df

df_benchmark = pd.concat(pd.read_parquet(i) for i in tqdm(snakemake.input['benchmark']))
df_benchmark = median_n_filter(df_benchmark)

if snakemake.params['subset_samples'] == True:
    df_benchmark['sample'] = df_benchmark['sample'].astype(str)
    df_benchmark = df_benchmark[df_benchmark['sample'].str.contains('CASE')]

# # rank methods
df_benchmark['outlier'] = df_benchmark['outlier'].fillna(False).astype(int)

models = [x for x in snakemake.params['model_dict'].values()]

for model in models:
    df_benchmark[model] = df_benchmark[model].fillna(0)
    df_benchmark[f'{model}_rank_mean'] = np.abs(df_benchmark[model].fillna(0)).rank(method='average', ascending=False)
    
for model in models:
    print(model)
    df_benchmark[f'{model}_tp_sum'] = df_benchmark.sort_values(by=model, ascending=False, key=abs)['outlier'].fillna(0).cumsum()

rank_cols = [f'{x}_rank_mean' for x in models]
tp_cols = [f'{x}_tp_sum' for x in models]

df_tp_rank_plot = df_benchmark[[*models, *rank_cols, *tp_cols]]

df_list = list()
for model in models:
    df_temp = df_tp_rank_plot[[model, f'{model}_rank_mean', f'{model}_tp_sum']]
    df_temp = df_temp.rename(columns = {model: 'model_score', f'{model}_rank_mean': 'rank_mean', f'{model}_tp_sum': 'tp_sum'})
    df_temp['model'] = model
    df_list.append(df_temp)
df_list = pd.concat(df_list)

# get points at defined cutoffs
model_cutoffs = snakemake.params['model_cutoffs']

df_defined_cutoffs = list()
for model in models:
    for cutoff, cutoff_value in model_cutoffs[model].items():
        df_defined_cutoffs.append(pd.DataFrame({
            'model': model,
            'cutoff': cutoff,
            'cutoff_value': cutoff_value,
            'tp_sum': df_benchmark[df_benchmark[model] >= cutoff_value].sort_values(by=model, ascending=True, key=abs).iloc[0][f'{model}_tp_sum'],
            'rank_mean': df_benchmark[df_benchmark[model] >= cutoff_value].sort_values(by=model, ascending=True, key=abs).iloc[0][f'{model}_rank_mean'],
        }, index=[0]))
df_defined_cutoffs = pd.concat(df_defined_cutoffs)

df_list.to_csv(snakemake.output['tp_among_all_preds_line'])
df_defined_cutoffs.to_csv(snakemake.output['tp_among_all_preds_discrete_cutoff'], index=False)

# # get means
df_list_mean = df_list.groupby(['model', 'tp_sum']).agg('mean').reset_index()
df_list_mean.to_csv(snakemake.output['tp_among_all_preds_line'].replace('.csv', '_interpolated.csv'), index=False)
