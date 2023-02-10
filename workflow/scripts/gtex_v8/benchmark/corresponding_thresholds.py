import pandas as pd
from tqdm import tqdm
import pickle
import os
import re
import numpy as np

performance = pickle.load(open(snakemake.input['aps_performance'], 'rb'))

# get ridd of last entry
df_list = list()
for model in performance.keys():
    df_list.append(
        pd.DataFrame({
            'precision': performance[model]['precision'][:-1],
            'recall': performance[model]['recall'][:-1],
            'threshold': performance[model]['threshold'],
            'model': model
        })
    )
df_pr = pd.concat(df_list)
df_pr.to_csv(snakemake.output['thresholds_per_model'], index=False)

# get corresponding thresholds
stepsize = 0.001
recall_cutoffs = np.arange(0+stepsize, 1, stepsize)
df_sel = list()

for recall in tqdm(recall_cutoffs):
    for model in performance.keys():
        _df = df_pr[
            (df_pr['model'] == model)
            & (df_pr['recall'] >= recall)
        ].sort_values(by='recall')
        df_sel.append({
            'model': model,
            'recall_rounded': recall,
            'recall': _df.iloc[0]['recall'],
            'precision': _df.iloc[0]['precision'],
            'threshold': _df.iloc[0]['threshold'],   
            }
        )
df_sel = pd.DataFrame(df_sel)
df_sel.to_csv(snakemake.output['thresholds'], index=False)
