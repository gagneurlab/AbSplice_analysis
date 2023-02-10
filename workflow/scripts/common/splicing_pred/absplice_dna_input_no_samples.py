import pandas as pd
from absplice import SplicingOutlierResult


# result = SplicingOutlierResult(
#     df_mmsplice = snakemake.input['pred_mmsplice'],
#     df_spliceai = snakemake.input['pred_spliceai'],
#     gene_map = snakemake.input['gene_map'],
#     gene_tpm = snakemake.input['gene_tpm'],
# )

# # result.absplice_dna_input.to_csv(snakemake.output['absplice_dna_input'])
# result.absplice_dna_input.to_parquet(snakemake.output['absplice_dna_input'], engine='pyarrow')

df_mmsplice = pd.read_parquet(snakemake.input['pred_mmsplice'], engine='pyarrow')
df_spliceai = pd.read_parquet(snakemake.input['pred_spliceai'], engine='pyarrow')

print(df_mmsplice.shape)
print(df_spliceai.shape)

if df_mmsplice.shape[0] == 0 or df_spliceai.shape[0] == 0:
    if df_mmsplice.shape[0] == 0:
        df_mmsplice = None

        from splicemap.splice_map import SpliceMap
        from tqdm import tqdm
        df5 = pd.concat(
            [SpliceMap.read_csv(i).df.assign(tissue=t) for i,t in tqdm(zip(snakemake.input['splicemap_5'], snakemake.params['gtex_tissues']))]
        )
        df3 = pd.concat(
            [SpliceMap.read_csv(i).df.assign(tissue=t) for i,t in tqdm(zip(snakemake.input['splicemap_3'], snakemake.params['gtex_tissues']))]
        )
        assert snakemake.wildcards['vcf_id'] not in df5['gene_id'].unique()
        assert snakemake.wildcards['vcf_id'] not in df3['gene_id'].unique()
    
    if df_spliceai.shape[0] == 0:
        df_spliceai = None
    
    result = SplicingOutlierResult(
        df_mmsplice = df_mmsplice,
        df_spliceai = df_spliceai,
        gene_map = snakemake.input['gene_map'],
        gene_tpm = snakemake.input['gene_tpm'],
    )
    
else:
    result = SplicingOutlierResult(
        df_mmsplice = snakemake.input['pred_mmsplice'],
        df_spliceai = snakemake.input['pred_spliceai'],
        gene_map = snakemake.input['gene_map'],
        gene_tpm = snakemake.input['gene_tpm'],
    )

# result.absplice_dna_input.to_csv(snakemake.output['absplice_dna_input'])
result.absplice_dna_input.to_parquet(snakemake.output['absplice_dna_input'], engine='pyarrow')
