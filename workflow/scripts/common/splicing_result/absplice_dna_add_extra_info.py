import pandas as pd

# read in dfs
df_absplice = pd.read_parquet(snakemake.input['absplice']).reset_index()
df_spliceai = pd.read_parquet(snakemake.input['spliceai'])
df_mmsplice_splicemap = pd.read_parquet(snakemake.input['mmsplice_splicemap'])

# Join extra info to absplice
index = ['variant', 'gene_id', 'tissue']
if 'sample' in df_absplice.columns:
    index = [*index, 'sample']

cols_spliceai = [
    'acceptor_gain', 'acceptor_loss', 'donor_gain', 'donor_loss', 
    'acceptor_gain_position', 'acceptor_loss_positiin', 'donor_gain_position', 'donor_loss_position', 
]
cols_mmsplice_splicemap = [
    'junction', 'event_type', 'splice_site', 'events', 
    'ref_psi', 'median_n'
]

df_absplice = df_absplice.set_index(index)
df_spliceai = df_spliceai.set_index(index)
df_mmsplice_splicemap = df_mmsplice_splicemap.set_index(index)

df_absplice_extra_info = df_absplice.join(
    df_spliceai[cols_spliceai]).join(
        df_mmsplice_splicemap[cols_mmsplice_splicemap])\
            .reset_index()
            
# save
df_absplice_extra_info.to_parquet(snakemake.output['absplice_extra_info'], index=False)
