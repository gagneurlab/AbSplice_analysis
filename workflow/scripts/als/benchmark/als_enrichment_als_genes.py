import pandas as pd
import numpy as np
from tqdm import tqdm

def median_n_cutoff(_df):
    # median_n cutoff for SpliceMap results
    _df.loc[
        _df['median_n_mmsplice_splicemap'] < snakemake.params['median_n_cutoff'], 
        'delta_logit_psi_mmsplice_splicemap'] = 0
    _df.loc[
        _df['median_n_mmsplice_splicemap_ref_psi'] < snakemake.params['median_n_cutoff'], 
        'delta_psi_mmsplice_splicemap_ref_psi'] = 0
    
    _df.loc[
        (_df['psi3_median_n_spliceai_splicemap'] < snakemake.params['median_n_cutoff']) \
        | (_df['psi5_median_n_spliceai_splicemap'] < snakemake.params['median_n_cutoff']),
        'delta_score_spliceai_splicemap'] = 0
    _df.loc[
        (_df['psi3_median_n_spliceai_splicemap_ref_psi'] < snakemake.params['median_n_cutoff']) \
        | (_df['psi5_median_n_spliceai_splicemap_ref_psi'] < snakemake.params['median_n_cutoff']),
        'delta_score_spliceai_splicemap_ref_psi'] = 0
    return _df

# read in benchmark
df_benchmark = pd.concat(
    [pd.read_parquet(i) for i in tqdm(snakemake.input['benchmark'])]
)

df_benchmark = df_benchmark.set_index(['gene_id','sample'])
df_benchmark = median_n_cutoff(df_benchmark)
df_benchmark['CASE'] = df_benchmark.index.get_level_values('sample').str.startswith('CASE')
df_benchmark = df_benchmark[df_benchmark['CASE'] == True]

# annotate als genes
df_als = pd.read_csv(snakemake.input['als_genes'])

df_gene_mapping = pd.read_csv('../../data/resources/common/hg38/gene_id_to_name_mapping_updated_with_gencode_gtf.tsv', sep='\t')
gene_map = dict(zip(df_gene_mapping['gene_name'], df_gene_mapping['gene_id']))
df_als['gene_id'] = df_als['gene_name'].map(gene_map)
df_als = df_als[~df_als['gene_id'].isna()]
als_genes = set(df_als['gene_id'])
num_als_genes = len(als_genes)

df_benchmark['als_gene'] = df_benchmark.index.get_level_values('gene_id').isin(als_genes)

# # Enrichment in als genes
model_cutoffs = snakemake.params['model_cutoffs']
cutoff_level = 'high'
models = [x for x in model_cutoffs.keys() if x in df_benchmark.columns]

df_category = list()
for model in tqdm(models):
    cutoff = model_cutoffs[model][cutoff_level]
    df_category.append(pd.DataFrame({
            'cutoff': cutoff,
            'model': model,
            'above_cutoff_all': df_benchmark[
                (np.abs(df_benchmark[model]) >= cutoff)].shape[0],
            'above_cutoff_in_category': df_benchmark[
                (np.abs(df_benchmark[model]) >= cutoff)
                & (df_benchmark['als_gene'] == True)].shape[0]
        }, index=[0]))
df_category = pd.concat(df_category)

# for random, get all protein coding genes
df_coding_genes = pd.read_csv(snakemake.input['coding_genes'])
num_all_coding_genes = df_coding_genes[df_coding_genes['gene_type'] == 'protein_coding'].shape[0]
df_random = pd.DataFrame.from_dict({
    'gene_category': ['als_gene'],
    'model': ['random'],
    'above_cutoff_all': [num_all_coding_genes],
    'above_cutoff_in_category': [num_als_genes]
})

df_enrichment = pd.concat([df_category, df_random])

df_enrichment['fraction_correct'] = df_enrichment['above_cutoff_in_category'] / df_enrichment['above_cutoff_all']
df_enrichment['enrichment'] = df_enrichment['fraction_correct'] / df_enrichment[df_enrichment['model'] == 'random']['fraction_correct'].values[0]

df_enrichment.to_csv(snakemake.output['enrichment_als_genes'], index=False)