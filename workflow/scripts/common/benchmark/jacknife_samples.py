import pandas as pd
from tqdm import tqdm
import numpy as np
from absplice_scripts.visualization.benchmark import jackknife_performance

def median_n_cutoff(df):
    # median_n cutoff for SpliceMap results
    df.loc[
        df['median_n_mmsplice_splicemap'] < snakemake.params['median_n_cutoff'], 
        'delta_logit_psi_mmsplice_splicemap'] = 0
    df.loc[
        df['median_n_mmsplice_splicemap'] < snakemake.params['median_n_cutoff'], 
        'delta_psi_mmsplice_splicemap'] = 0

    df.loc[
        (df['psi3_median_n_spliceai_splicemap'] < snakemake.params['median_n_cutoff']) | (df['psi5_median_n_spliceai_splicemap'] < snakemake.params['median_n_cutoff']),
        'delta_score_spliceai_splicemap'] = 0
    df.loc[
        (df['psi3_median_n_spliceai_splicemap_ref_psi'] < snakemake.params['median_n_cutoff']) | (df['psi5_median_n_spliceai_splicemap_ref_psi'] < snakemake.params['median_n_cutoff']),
        'delta_score_spliceai_splicemap_ref_psi'] = 0
    
    return df


df_benchmark = pd.concat(
    [pd.read_parquet(i) for i in tqdm(snakemake.input['df_benchmark'])]
)

if 'samples' in snakemake.params.keys():
    df_benchmark = df_benchmark[df_benchmark['sample'].isin(snakemake.params['samples'])]
df_benchmark = df_benchmark.set_index(['gene_id', 'sample'])
df_benchmark = median_n_cutoff(df_benchmark)

df_auPRC, df_stats = jackknife_performance(
    df_benchmark, 
    model_dict = snakemake.params['model_dict'], 
    outlier_column = snakemake.params['outlier_column']
)

df_auPRC.to_csv(snakemake.output['auPRC'], index=False)
df_stats.to_csv(snakemake.output['stats'], index=False)