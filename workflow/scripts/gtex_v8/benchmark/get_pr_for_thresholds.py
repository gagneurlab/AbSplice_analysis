import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import os
from pathlib import Path
import pickle
from tqdm import tqdm

# def match_model_cutoff_to_recall(df_selected_cutoffs, model_compare='SpliceAI', model='AbSplice_DNA', threshold=0.5):
#     recall_rounded = df_selected_cutoffs[
#         (df_selected_cutoffs['model'] == model_compare)
#         & (df_selected_cutoffs['threshold'] >= threshold)
#     ].sort_values(by='threshold').iloc[0]['recall_rounded']
    
#     model_cutoff = df_selected_cutoffs[
#         (df_selected_cutoffs['model'] == model)
#         & (df_selected_cutoffs['recall_rounded'] == recall_rounded)
#     ]['threshold'].values[0]
    
#     return model_cutoff

model_cutoffs = {
    'SpliceAI': {
        'high': 0.8,
        'medium': 0.5,
        'low': 0.2,
    },
    'MMSplice': {
        'high': 2,
        'medium': 1.5,
        'low': 1,
    },
    'MMSplice SpliceMap + Ref_PSI': {
        'high': 0.2,
        'medium': 0.1,
        'low': 0.05,
    },
    'AbSplice-DNA': {
        'high': 0.2,
       'medium': 0.05,
       'low': 0.01
    }
}
model_cutoffs['SpliceAI SpliceMap'] = model_cutoffs['SpliceAI']
model_cutoffs['SpliceAI SpliceMap + Ref PSI'] = model_cutoffs['SpliceAI']
model_cutoffs['MMSplice SpliceMap'] = model_cutoffs['MMSplice']
models = model_cutoffs.keys()

df_pr = pd.read_csv(snakemake.input['thresholds_per_model'])

def get_pr_for_threshold(df_pr, threshold, model):
    _df = df_pr[df_pr['model'] == model]
    _df = _df.sort_values(by='threshold')
    _df = _df[_df['threshold'] >= threshold]
    precision = _df.iloc[0]['precision']
    recall = _df.iloc[0]['recall']
    return precision, recall

df_defined_cutoffs = list()
for model in models:
    if model in set(df_pr['model']):
        for cutoff, cutoff_value in model_cutoffs[model].items():
            precision, recall = get_pr_for_threshold(df_pr, cutoff_value, model)
            df_defined_cutoffs.append(pd.DataFrame({
                'model': model,
                'cutoff': cutoff,
                'cutoff_value': cutoff_value,
                'precision': precision,
                'recall': recall,
            }, index=[0]))
df_defined_cutoffs = pd.concat(df_defined_cutoffs)
df_defined_cutoffs.to_csv(snakemake.output['threshold_points_pr_curve'], index=False)