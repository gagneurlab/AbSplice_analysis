import pandas as pd
import pyranges as pr
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors import variants_to_pyranges

def cleanup_df(df, columns):
    df = df[columns]
    df['gene_id'] = df['gene_id'].apply(lambda x: x.split('.')[0])
    df['variant'] = df['variant'].apply(lambda x: x. __str__())
    df = df.set_index(['variant', 'gene_id'])
    return df

# Variants
df_vars = pd.read_csv(snakemake.input['variants'])
v = Variant
v_list = df_vars['variant'].apply(lambda x: v.from_str(x)).values
v_list = list(set(v_list))
pr_vars = variants_to_pyranges(v_list)

# genes and exons from GTF
pr_gtf_complete = pr.read_gtf(snakemake.input['gtf'])
if 'gene_type' not in pr_gtf_complete.columns and 'gene_biotype' in pr_gtf_complete.columns:
    pr_gtf_complete = pr.PyRanges(pr_gtf_complete.df.rename(columns={'gene_biotype': 'gene_type'}))
pr_gtf_exons = pr_gtf_complete[(pr_gtf_complete.Feature == 'exon') & (pr_gtf_complete.gene_type == 'protein_coding')]
pr_gtf_genes = pr_gtf_complete[(pr_gtf_complete.Feature == 'gene') & (pr_gtf_complete.gene_type == 'protein_coding')]

# get distance of variants to nearest exon
pr_vars_exon_dist = pr_vars.k_nearest(pr_gtf_exons, k=1, overlap=True, ties='first') # exonic distance
df_vars_exon_dist = cleanup_df(pr_vars_exon_dist.df, columns = ['variant', 'gene_id', 'exon_id', 'exon_number', 'Distance']) # clean up df
df_vars_exon_dist['exonic'] = df_vars_exon_dist['Distance'] == 0

pr_vars_intron_dist = pr_vars.k_nearest(pr_gtf_exons, k=1, overlap=False, ties='first') # intronic distance
df_vars_intron_dist = cleanup_df(pr_vars_intron_dist.df, columns = ['variant', 'gene_id', 'exon_id', 'exon_number', 'Distance']) # clean up df

df_vars_intron_exon_dist = df_vars_intron_dist.join(df_vars_exon_dist[['exonic', 'Distance']], rsuffix='_exonic')

# get gene_id for variants
pr_vars_gene = pr_vars.join(pr_gtf_genes) #default is None, i.e. inner join
df_vars_gene = cleanup_df(pr_vars_gene.df, columns=['variant', 'gene_id', 'gene_name', 'gene_type']) # clean up df

# only consider variants within gene boundaries (i.e. inner join)
df_vars_annotated = df_vars_intron_exon_dist.join(df_vars_gene, how='inner', rsuffix='_gene')
df_vars_annotated = df_vars_annotated.reset_index()

# recover rest of annotations from variants
df_vars_annotated = df_vars_annotated.set_index('variant').join(df_vars.set_index('variant'), how='left').reset_index()
df_vars_annotated.to_csv(snakemake.output['variants_annotated'], index=False)