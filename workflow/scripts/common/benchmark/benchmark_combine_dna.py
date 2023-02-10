import pandas as pd
from absplice.utils import get_abs_max_rows

params_dict = dict()
for k,v in snakemake.params.items():
    if 'col' in k:
        params_dict[k] = v
        
# universe (gene, sample) with at least one rare variant
universe = pd.read_csv(snakemake.input['universe']).set_index(['gene_id', 'sample'])
universe['tissue'] = snakemake.params['tissue']

# outliers
outliers = pd.read_csv(snakemake.input['outliers'])
outliers['outlier'] = True
outliers = outliers[['gene_id', 'sample', 'outlier']]
outliers = outliers.drop_duplicates().set_index(['gene_id', 'sample']) # there can be multiple outliers on a gene. drop duplicates

# create benchmark universe
df_benchmark = universe.join(outliers, how='left')

# predictions
pred_tools = snakemake.params['pred_tools']
tissue_specific_tools = snakemake.params['tissue_specific_tools']

for pred_tool in pred_tools:
    try:
        df = pd.read_parquet(snakemake.input[f'{pred_tool}'])
    except:
        df = pd.read_csv(snakemake.input[f'{pred_tool}'])
    df = df.reset_index()
    df = df\
        .set_index(['gene_id', 'sample'])[params_dict[f'cols_{pred_tool}']]\
        .add_suffix(f'_{pred_tool}')
    if 'tissue_subset' not in snakemake.params.keys(): # TODO: This is a quick fix for tissue subset. 
        if pred_tool in tissue_specific_tools:
            if snakemake.params['splicemaps'] == 'gtex':
                df = df[df[f'tissue_{pred_tool}'] == snakemake.params['gtex_tissue']]
            elif snakemake.params['splicemaps'] == 'dataset':
                df = df[df[f'tissue_{pred_tool}'] == snakemake.params['tissue']]
    df_benchmark = df_benchmark.join(df, how='left')
    
# save benchmark
df_benchmark = df_benchmark.reset_index()
df_benchmark.to_parquet(snakemake.output['combined_benchmark'], index=False)