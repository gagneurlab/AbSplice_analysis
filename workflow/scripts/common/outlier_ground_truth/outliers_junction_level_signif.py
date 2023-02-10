import pandas as pd
from kipoiseq.extractors.vcf_seq import MultiSampleVCF
from splicemap.count_table import SpliceCountTable as CountTable
from absplice_scripts.data.DROP_annotations import sample_individual_mapping, annotate_junctions

# Initialization (all contain 'chr')
ct = CountTable.read_csv(snakemake.input['count_table_updated'])

df_junction_level = pd.read_csv(snakemake.input['junction_level'], sep='\t')\
                    .rename(columns={'sampleID': 'sample'})
df_junction_level = annotate_junctions(df_junction_level)

df_anno = pd.read_csv(snakemake.input['sample_map'], sep='\t')
samples = set(df_anno[
    (df_anno['DROP_GROUP'].str.contains(snakemake.params['tissue_DROP']))
    & (~df_anno[snakemake.params['value_assay']].isna())
    & (~df_anno[snakemake.params['key_assay']].isna())
][snakemake.params['value_assay']])

# Filtering
sample_map = sample_individual_mapping(
    annotation_table_path=snakemake.input['sample_map'],
    key_assay=snakemake.params['key_assay'],                           
    value_assay=snakemake.params['value_assay'])

# filter for samples with WGS
df_junction_level['sample'] = df_junction_level['sample'].astype(str)
df_junction_level['sample'] = df_junction_level['sample'].map(sample_map)
df_junction_level = df_junction_level[df_junction_level['sample'].isin(samples)]

# filter for protein coding
    # df_ct_annotation = pd.read_csv(snakemake.input['count_table_with_annotation'])
    # df_junction_level = df_junction_level.set_index('junctions')\
    #                 .join(df_ct_annotation.set_index('junctions')[['gene_id', 'gene_name', 'gene_type']]).reset_index()
    # df_junction_level = df_junction_level[df_junction_level['gene_type'] == 'protein_coding']
gene_map = pd.read_csv(snakemake.input['gene_map'], sep='\t')
gene_map = dict(zip(gene_map['gene_name'], gene_map['gene_id']))
df_junction_level['gene_id'] = df_junction_level['hgncSymbol'].map(gene_map)

df_coding = pd.read_csv(snakemake.input['coding_genes'])
df_junction_level = df_junction_level[
    df_junction_level['gene_id'].isin(set(df_coding['gene_id']))
]


# filter by cutoffs
df_junction_level = df_junction_level[
    (df_junction_level['padjust'] <= snakemake.params['padjust_junction_cutoff'])
    & (df_junction_level['totalCounts'] >= snakemake.params['totalCounts_cutoff'])
#     & (df_junction_level['deltaPsi'].abs() >= snakemake.params['delta_psi_cutoff'])
]

df_junction_psi5 = df_junction_level[df_junction_level['type'] == 'psi5'].set_index('junctions')
df_junction_psi3 = df_junction_level[df_junction_level['type'] == 'psi3'].set_index('junctions')
df_junction_theta = df_junction_level[df_junction_level['type'] == 'theta']

df_junction_psi5 = df_junction_psi5.join(ct.splice_site5).join(ct.event5)
df_junction_psi3 = df_junction_psi3.join(ct.splice_site3).join(ct.event3)

df_junction_theta['events'] = df_junction_theta['junctions']
df_junction_theta = df_junction_theta.set_index('junctions')

pd.concat([
    df_junction_psi5, 
    df_junction_psi3,
    df_junction_theta
]).to_csv(snakemake.output['junction_signif'])
