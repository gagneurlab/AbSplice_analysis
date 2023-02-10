import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import os
from pathlib import Path
import pickle
from tqdm import tqdm

from interpret.glassbox import ExplainableBoostingClassifier
from interpret.provider import InlineProvider
from interpret import set_visualize_provider
from interpret import show, show_link
from interpret import preserve
from interpret import get_show_addr

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GroupKFold

model_params = {
    "max_depth": 3,
    "n_estimators": 100,
    "criterion": "entropy",
    "random_state": 43,
    "class_weight": "balanced",
    "max_features": "auto"
}


df_absplice_input = pd.read_csv(snakemake.input['absplice_input'])
df_absplice_input['splice_site_is_expressed'] = (df_absplice_input['median_n'] > snakemake.params['median_n_cutoff']).astype(int)
df_absplice_input['gene_is_expressed'] = (df_absplice_input['gene_tpm'] > snakemake.params['gene_tpm_cutoff']).astype(int)

features = snakemake.params['feature_string'].split('__')

df_ensemble_all_gtex = df_absplice_input.fillna(0)

if snakemake.params['abs_features'] == 'True':
    abs_feature = True
elif snakemake.params['abs_features'] == 'False':
    abs_feature = False
else:
    raise ValueError('abs_features needs to be boolean')
if abs_feature:
    df_ensemble_all_gtex[features] = np.abs(df_ensemble_all_gtex[features])
    
df_ensemble_all_gtex['outlier'] = df_ensemble_all_gtex['outlier'].astype(int)

X_train = df_ensemble_all_gtex[features]
y_train = df_ensemble_all_gtex['outlier']

if snakemake.params['classifier'] == 'interpretml':
    model_all_gtex = ExplainableBoostingClassifier()
elif snakemake.params['classifier'] == 'califorest':
    model_all_gtex = CaliForest()
    X_train = X_train.values
    y_train = y_train.values
elif snakemake.params['classifier'] == 'random_forest':
    X_train = X_train.values
    y_train = y_train.values
    model_all_gtex = RandomForestClassifier(**model_params)
elif snakemake.params['classifier'] == 'logistic_regression':
    X_train = X_train.values
    y_train = y_train.values
    model_all_gtex = LogisticRegression(random_state=0)
model_all_gtex.fit(X_train, y_train)

save_dir = os.path.dirname(snakemake.output['absplice_whole_GTEx'])
path = Path(save_dir)
path.mkdir(parents=True, exist_ok=True)

print(snakemake.output['absplice_whole_GTEx'])
with open(snakemake.output['absplice_whole_GTEx'], 'wb') as file:
    pickle.dump(model_all_gtex, file)


