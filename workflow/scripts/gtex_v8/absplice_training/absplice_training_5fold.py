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

def train_model_classification(
    df, save_dir, 
    features,
    chosen_model='interpretml',
    model_params=model_params,
    abs_features=True):
    
    df = df.fillna(0)
    df['outlier'] = df['outlier'].astype(int)
    
    X = df.set_index(['variant', 'gene_id', 'sample', 'tissue'])[features]
    y = df.set_index(['variant', 'gene_id', 'sample', 'tissue'])['outlier']
    print(f'y shape: {y.shape}')
    print(f'y outliers: {y.sum()}')

    groups = df['sample']
    gkf = GroupKFold(n_splits=5)
    results_all = list()
    models_all = list()
    
    fold = 0
    for train, test in gkf.split(X, y, groups=groups):

        if chosen_model == 'interpretml':
            print('interpretml used')
            X_train = X[features].iloc[train]
            X_test = X[features].iloc[test]
            y_train = y.iloc[train]
            y_test = y.iloc[test]
            if abs_features:
                X_train = np.abs(X_train)
                X_test = np.abs(X_test)
            model = ExplainableBoostingClassifier()
            
        elif chosen_model == 'califorest':
            print('califorest used')
            X_train = X[features].iloc[train].values
            X_test = X[features].iloc[test].values
            y_train = y.iloc[train].values
            y_test = y.iloc[test].values
            if abs_features:
                X_train = np.abs(X_train)
                X_test = np.abs(X_test)
            model = CaliForest()
            
        elif chosen_model == 'random_forest':
            print('random forest used')
            X_train = X[features].iloc[train].values
            X_test = X[features].iloc[test].values
            y_train = y.iloc[train].values
            y_test = y.iloc[test].values
            if abs_features:
                X_train = np.abs(X_train)
                X_test = np.abs(X_test)
            model = RandomForestClassifier(**model_params)
            
        elif chosen_model == 'logistic_regression':
            print('logistic regression used')
            X_train = X[features].iloc[train].values
            X_test = X[features].iloc[test].values
            y_train = y.iloc[train].values
            y_test = y.iloc[test].values
            if abs_features:
                X_train = np.abs(X_train)
                X_test = np.abs(X_test)
            model = LogisticRegression(random_state=0)

        model.fit(X_train, y_train)
        models_all.append(model)
        y_pred = model.predict_proba(X_test)[:, 1]
        
        X_test_all = np.abs(X[features].iloc[test])
        
        pickle_filename = os.path.join(save_dir, 'cross_val=' + str(fold) + '.pkl')
        with open(pickle_filename, 'wb') as file:
            pickle.dump(model, file)

        results = pd.DataFrame({
            'variant': y.iloc[test].index.get_level_values('variant').values,
            'gene_id': y.iloc[test].index.get_level_values('gene_id').values,
            'sample': y.iloc[test].index.get_level_values('sample').values,
            'tissue': y.iloc[test].index.get_level_values('tissue').values,
            'y_pred': y_pred,
            'y_test': y_test,
            'fold': fold})

        for feature in features:
            results[feature] = X_test_all.iloc[:, features.index(feature)].values
            
        results_all.append(results)
        fold += 1

        del model

    results_all_df = pd.concat(results_all) 
    
    results_filename = os.path.join(save_dir, 'results_all.csv')
    results_all_df.to_csv(results_filename, index=False)
    
    return results_all_df, models_all


df_absplice_input = pd.read_csv(snakemake.input['absplice_input'])
df_absplice_input['splice_site_is_expressed'] = (df_absplice_input['median_n'] > snakemake.params['median_n_cutoff']).astype(int)
df_absplice_input['gene_is_expressed'] = (df_absplice_input['gene_tpm'] > snakemake.params['gene_tpm_cutoff']).astype(int)

features = snakemake.params['feature_string'].split('__')

save_dir = os.path.dirname(snakemake.output['absplice_5fold_crossval'])
path = Path(save_dir)
path.mkdir(parents=True, exist_ok=True)

if snakemake.params['abs_features'] == 'True':
    abs_feature = True
elif snakemake.params['abs_features'] == 'False':
    abs_feature = False
else:
    raise ValueError('abs_features needs to be boolean')

results_all, model_all = train_model_classification(
    df=df_absplice_input, 
    save_dir = save_dir,
    features=features, 
    chosen_model=snakemake.params['classifier'],
    model_params=model_params,
    abs_features=abs_feature
)


