import pandas as pd
import pickle
from tqdm import tqdm
import re

wildcards = list()
results = list()

for path in tqdm(snakemake.input['benchmark']):
    wildcard = re.findall('^.+performance_tissue=(.+)_aps.pkl', path)[0]
    benchmark = pickle.load(open(path, 'rb'))
    for model, result in benchmark.items():
        results.append((model, result['aps']))
        wildcards.append(wildcard)

df_benchmark = pd.concat([
    pd.DataFrame(results, columns=['model', 'Average Precision Score']),
    pd.DataFrame(wildcards, columns=['tissue'])
], axis=1)

scores_mean = pd.DataFrame(df_benchmark.groupby('model').mean())['Average Precision Score'].to_dict()
scores_std = pd.DataFrame(df_benchmark.groupby('model').std())['Average Precision Score'].to_dict()

df_benchmark['mean'] = df_benchmark['model'].map(lambda x: scores_mean[x])
df_benchmark['std'] = df_benchmark['model'].map(lambda x: scores_std[x])

df_benchmark.to_csv(snakemake.output['df_performance_boxplot'], index=False)
