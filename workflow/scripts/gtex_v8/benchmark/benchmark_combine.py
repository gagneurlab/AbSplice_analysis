import pandas as pd
from tqdm import tqdm
from absplice.utils import get_abs_max_rows

index = snakemake.params['unique_index']

# # update cols for absplice
# features = snakemake.wildcards['feature_string'].split('__')
params_dict = dict()
for k,v in snakemake.params.items():
    if 'col' in k:
        params_dict[k] = v
# params_dict['cols_absplice'] = [*params_dict['cols_absplice'], *features]

# universe (gene, sample) with at least one rare variant
universe = pd.read_csv(snakemake.input['universe']).set_index(index)
universe['tissue'] = snakemake.wildcards['tissue']

# outliers
outliers = pd.read_csv(snakemake.input['outliers']).set_index(index)
outliers['outlier'] = True

# create benchmark universe
df_benchmark = universe.join(outliers, how='left')

# predictions
pred_tools = snakemake.params['pred_tools']
tissue_specific_tools = snakemake.params['tissue_specific_tools']

for pred_tool in tqdm(pred_tools):
    if pred_tool in snakemake.input.keys():
        try:
            df = pd.read_parquet(snakemake.input[f'{pred_tool}'])
        except:
            df = pd.read_csv(snakemake.input[f'{pred_tool}'])
        df = df.reset_index()
        df = df\
            .set_index(index)[params_dict[f'cols_{pred_tool}']]\
            .add_suffix(f'_{pred_tool}')
        if pred_tool in tissue_specific_tools:
            df = df[df[f'tissue_{pred_tool}'] == snakemake.wildcards['tissue']]
        df_benchmark = df_benchmark.join(df, how='left')

# save benchmark
df_benchmark = df_benchmark.reset_index()
df_benchmark.to_parquet(snakemake.output['combined_benchmark'], index=False)
