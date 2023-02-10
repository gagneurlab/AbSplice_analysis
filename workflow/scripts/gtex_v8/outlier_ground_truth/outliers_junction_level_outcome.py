import pandas as pd
from absplice_scripts.data.DROP_annotations import sample_individual_mapping, annotate_junctions

df_ct_annotation = pd.read_csv(snakemake.input['count_table_with_annotation'])
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
df_junction_level = df_junction_level.set_index('junctions')\
                .join(df_ct_annotation.set_index('junctions')[['gene_id', 'gene_name', 'gene_type']]).reset_index()
df_junction_level = df_junction_level[df_junction_level['gene_type'] == 'protein_coding']

# theta outliers (not in outcome annotation)
df_outliers_signif = pd.read_csv(snakemake.input['outliers_signif'])

# Psi outliers
df_psi = df_outliers_signif[
    (df_outliers_signif['type_j'] == 'psi5')
    | (df_outliers_signif['type_j'] == 'psi3')
]

if df_psi.shape[0] > 0:
    df_psi.loc[df_psi['deltaPsi_j'] >= 0, 'aberrantSpliceType'] = 'psi_increase'
    df_psi.loc[df_psi['deltaPsi_j'] < 0, 'aberrantSpliceType'] = 'psi_decrease'

cols = [x for x in df_psi.columns if '_g' not in x]
df_psi = df_psi[cols]

replace_dict = {i:i.replace('_j', '') for i in cols}
df_psi = df_psi.rename(columns=replace_dict)

# theta outliers
df_theta = df_outliers_signif[
    (df_outliers_signif['type_j'] == 'theta')
]

if df_theta.shape[0] > 0:
    df_theta.loc[df_theta['deltaPsi_j'] >= 0, 'aberrantSpliceType'] = 'theta_increase'
    df_theta.loc[df_theta['deltaPsi_j'] < 0, 'aberrantSpliceType'] = 'theta_decrease'

cols = [x for x in df_theta.columns if '_g' not in x]
df_theta = df_theta[cols]

replace_dict = {i:i.replace('_j', '') for i in cols}
df_theta = df_theta.rename(columns=replace_dict)

# join all together
df_junction_level = pd.concat([df_junction_level, df_theta, df_psi])
df_junction_level.to_csv(snakemake.output['outlier_outcome'], index=False)
