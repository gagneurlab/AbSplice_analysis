import pandas as pd
from kipoiseq.extractors.vcf_seq import MultiSampleVCF
from absplice_scripts.data.DROP_annotations import sample_individual_mapping, annotate_junctions

df_anno = pd.read_csv(snakemake.input['sample_map'], sep='\t')
print(snakemake.params['tissue_DROP'])
samples = set(df_anno[
    (df_anno['DROP_GROUP'].str.contains(snakemake.params['tissue_DROP']))
    & (~df_anno[snakemake.params['value_assay']].isna())
    & (~df_anno[snakemake.params['key_assay']].isna())
][snakemake.params['value_assay']])

# sample map
sample_map = sample_individual_mapping(
    annotation_table_path=snakemake.input['sample_map'],
    key_assay=snakemake.params['key_assay'],                           
    value_assay=snakemake.params['value_assay'])

# gene level outliers
df_gene_level = pd.read_csv(snakemake.input['gene_level'], sep='\t')\
                    .rename(columns={'sampleID': 'sample'})
df_gene_level = annotate_junctions(df_gene_level)

# filter for samples with WGS
df_gene_level['sample'] = df_gene_level['sample'].astype(str)
df_gene_level['sample'] = df_gene_level['sample'].map(sample_map)
df_gene_level = df_gene_level[df_gene_level['sample'].isin(samples)]

# filter for protein coding
    # # count table with inferred gene annotation
    # df_ct_annotation = pd.read_csv(snakemake.input['count_table_with_annotation'])
    # df_gene_level = df_gene_level.set_index('junctions')\
    #         .join(df_ct_annotation.set_index('junctions')[['gene_id', 'gene_name', 'gene_type']]).reset_index()
    # df_gene_level = df_gene_level[df_gene_level['gene_type'] == 'protein_coding']
gene_map = pd.read_csv(snakemake.input['gene_map'], sep='\t')
gene_map = dict(zip(gene_map['gene_name'], gene_map['gene_id']))
df_gene_level['gene_id'] = df_gene_level['hgncSymbol'].map(gene_map)

df_coding = pd.read_csv(snakemake.input['coding_genes'])
df_gene_level = df_gene_level[
    df_gene_level['gene_id'].isin(set(df_coding['gene_id']))
]

# filter by cutoffs
df_gene_level = df_gene_level[
    df_gene_level['padjustGene'] <= snakemake.params['padjustGene_cutoff']
]

# save results
assert df_gene_level.shape[0] > 0
df_gene_level.to_csv(snakemake.output['result'], index=False)
