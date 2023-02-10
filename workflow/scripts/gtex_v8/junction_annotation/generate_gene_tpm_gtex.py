import pandas as pd

df = pd.read_csv(snakemake.input['gene_expression_raw'], skiprows=2, sep='\t')

df = df.rename(columns={
    'Name': 'gene_id', 
    'Description': 'gene_name'})

df['gene_id_orig'] = df['gene_id']
df['PAR_Y'] = df['gene_id'].apply(lambda x: 'PAR_Y' in x)
df = df[df['PAR_Y'] == False]
df['gene_id'] = df['gene_id'].apply(lambda x: x.split('.')[0])

df = df.rename(columns=snakemake.params['tissue_map'])
df.to_csv(snakemake.output['gene_tpm_wide'], index=False)

df = df.drop(columns=['gene_name', 'gene_id_orig', 'PAR_Y'])
df = df.melt(id_vars='gene_id', var_name='tissue', value_name='gene_tpm')
df.to_csv(snakemake.output['gene_tpm'], index=False)



