import pandas as pd
import pyranges as pr


# calculate TPM
gr = pr.read_gtf(snakemake.input['gtf_file'])
df_gene = gr[gr.Feature == 'gene'].df
df_gene['length'] = df_gene['End'] - df_gene['Start']
df_gene['gene_id'] = df_gene['gene_id'].str.split('.').str.get(0)
df_gene = df_gene.set_index('gene_id')[['length']]

df_gene_count = pd.read_csv(snakemake.input['gene_count'])
df_gene_count = df_gene_count.rename(columns={'Unnamed: 0': 'gene_id'}).set_index('gene_id')
samples = df_gene_count.columns[df_gene_count.columns.str.startswith('CTRL')]
df_gene_count = df_gene_count[samples]

df_gene_count = df_gene_count.join(df_gene, how='inner')
df_gene_count = df_gene_count.div(df_gene_count['length'], axis=0)
del df_gene_count['length']
df_gene_count = df_gene_count.div(df_gene_count.sum(axis=0) / 1_000_000, axis=1)
df_gene_count.to_csv(snakemake.output['gene_count'], index=False)

# get median values of expression
df_gene_tpm = df_gene_count.median(axis=1)
df_gene_tpm = pd.DataFrame(df_gene_tpm, columns=[snakemake.params['tissue']])

df_gene_map = gr[gr.Feature == 'gene'].df[['gene_id', 'gene_name']]
df_gene_map['gene_id'] = df_gene_map['gene_id'].apply(lambda x: x.split('.')[0])
gene_map = dict(zip(df_gene_map['gene_id'], df_gene_map['gene_name']))

df_gene_tpm = df_gene_tpm.reset_index()
df_gene_tpm['gene_name'] = df_gene_tpm['gene_id'].map(gene_map)
df_gene_tpm = df_gene_tpm[['gene_id', 'gene_name', snakemake.params['tissue']]]
df_gene_tpm = df_gene_tpm.drop_duplicates()
assert df_gene_tpm.shape[0] == len(set(df_gene_tpm['gene_id']))
df_gene_tpm.to_csv(snakemake.output['gene_tpm_wide'], index=False)

df_gene_tpm = df_gene_tpm[['gene_id', snakemake.params['tissue']]]
df_gene_tpm = df_gene_tpm.rename(columns={snakemake.params['tissue']: 'gene_tpm'})
df_gene_tpm['tissue'] = snakemake.params['tissue']
df_gene_tpm.to_csv(snakemake.output['gene_tpm'], index=False)

