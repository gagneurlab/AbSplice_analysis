import pandas as pd
from tqdm import tqdm
import os
import re
from absplice_scripts.visualization.benchmark import get_performance

def median_cutoff(df):
    # median_n cutoff for SpliceMap results
    df.loc[
        df['median_n_mmsplice_splicemap'] < snakemake.params['median_n_cutoff'], 
        'delta_logit_psi_mmsplice_splicemap'] = 0
    df.loc[
        df['median_n_mmsplice_splicemap_ref_psi'] < snakemake.params['median_n_cutoff'], 
        'delta_psi_mmsplice_splicemap_ref_psi'] = 0

    df.loc[
        (df['psi3_median_n_spliceai_splicemap'] < snakemake.params['median_n_cutoff']) \
        | (df['psi5_median_n_spliceai_splicemap'] < snakemake.params['median_n_cutoff']),
        'delta_score_spliceai_splicemap'] = 0
    df.loc[
        (df['psi3_median_n_spliceai_splicemap_ref_psi'] < snakemake.params['median_n_cutoff']) \
        | (df['psi5_median_n_spliceai_splicemap_ref_psi'] < snakemake.params['median_n_cutoff']),
        'delta_score_spliceai_splicemap_ref_psi'] = 0
    return df


def filter_var_category(df, var_category):
    if var_category != 'all':
        df = df[df[f'{var_category}'] == True]
    return df


def filter_outlier_category(df, outlier_category):
    if outlier_category != 'all':
        if df['outlier'].sum() > 0:
            if outlier_category == 'splicing_efficiency_change':
                df = df[
                    (
                        (df['outlier'] == True)
                        & (
                            (df['theta_decrease'] == True)
                            | (df['theta_increase'] == True)
                        )
                    )
                    | (
                        df['outlier'] != True
                    )
                ]
                 
            elif outlier_category == 'any_psi5_psi3_change':
                df = df[
                    (
                        (df['outlier'] == True)
                        & (
                            (df['psi_decrease'] == True)
                            | (df['psi_increase'] == True)
                        )
                    )
                    | (
                        df['outlier'] != True
                    )
                ]
                 
            else:
                df = df[
                    (df[f'{outlier_category}'] == True)
                    | (df['outlier'] != True)
                ] 
        else:
            print('File contains no outliers.')
            assert df['outlier'].sum() == 0
    return df



df_benchmark = list()
for i in tqdm(snakemake.input['benchmark']):
    df = pd.read_parquet(i)
    if df.shape[0] > 0:
        # filter for categories
        df = filter_var_category(df, var_category=snakemake.wildcards['var_category'])
        df = filter_outlier_category(df, outlier_category=snakemake.wildcards['outlier_category'])
        # filter for median_n
        df = median_cutoff(df)
        df_benchmark.append(df)
    else:
        print(f'DF IS EMPTY. \naffected file: {i}')
df_benchmark = pd.concat(df_benchmark)

# save performance
df_performance, performance = get_performance(
    df = df_benchmark, 
    outlier_column = snakemake.params['outlier_column'], 
    model_dict = snakemake.params['model_dict'],
)

df_performance.to_csv(snakemake.output['df_performance'], index=False)

# save statistics of outliers
df_stats = dict()
df_stats = pd.DataFrame({
    'var_category': snakemake.wildcards['var_category'],
    'outlier_category': snakemake.wildcards['outlier_category'],
    'num_outliers': df_benchmark[df_benchmark['outlier'] == True].shape[0],
    'num_non_outliers': df_benchmark[df_benchmark['outlier'] != True].shape[0],
    'num_benchmark_df': df_benchmark.shape[0],
}, index=[0])

df_stats.to_csv(snakemake.output['benchmark_stats'], index=False)