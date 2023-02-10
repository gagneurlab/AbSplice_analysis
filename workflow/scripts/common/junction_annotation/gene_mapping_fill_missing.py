import pandas as pd

df_gene_mapping = pd.read_csv(snakemake.input['gene_mapping_raw'], sep='\t')
df_gene_mapping = df_gene_mapping.rename(
    columns={'Ensembl_gene_ID': 'gene_id', 'any_symbol': 'gene_name'})


df_coding_genes = pd.read_csv(snakemake.input['coding_genes'])[['gene_id', 'gene_name']]
df_coding_genes['symbol_origin'] = 'Approved_symbol'

missing_genes = set(df_coding_genes['gene_id']).difference(set(df_gene_mapping['gene_id']))
df_missing_genes = df_coding_genes[df_coding_genes['gene_id'].isin(missing_genes)]

df_gene_mapping = pd.concat([
    df_gene_mapping, df_missing_genes
])

df_gene_mapping.to_csv(snakemake.output['gene_mapping'], sep='\t', index=False)