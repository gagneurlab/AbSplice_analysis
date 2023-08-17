import pandas as pd
import pyranges as pr
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors import variants_to_pyranges

# GTF
gr_gtf = pr.read_gtf(snakemake.input['gtf'])
# subset protein coding
gr_gtf = gr_gtf.subset(lambda df: df['gene_type'] == 'protein_coding')

# clean gene_ids in gtf
df_gtf = gr_gtf.df
df_gtf['gene_id'] = df_gtf['gene_id'].apply(lambda x: x.split('.')[0])
gr_gtf = pr.PyRanges(df_gtf)

# subset genes
pr_genes = gr_gtf.subset(lambda df: df['Feature'] == 'gene')

# SQUIRLS
df_squirls = pd.read_csv(snakemake.input['model'])

# get unique variants into pyranges
v_list = df_squirls['variant'].apply(lambda x: Variant.from_str(x)).values
v_list = list(set(v_list))
pr_vars = variants_to_pyranges(v_list)

# Get overlap
df_gene_var = pr_genes.join(pr_vars).df
df_gene_var = df_gene_var[['variant', 'gene_id']].drop_duplicates()
df_gene_var['variant'] = df_gene_var['variant'].astype(str)

df_squirls = df_squirls.set_index('variant').join(
    df_gene_var.set_index('variant')).reset_index()

df_squirls = df_squirls[['variant', 'gene_id', 'squirls_scores']].set_index('variant')

# annotate samples
if 'var_samples_df' in snakemake.input.keys():
    df_vcf_annotation = pd.read_csv(snakemake.input['var_samples_df'])\
                .set_index('variant')

    df_squirls = df_squirls.join(
        df_vcf_annotation, how='inner')\
        .reset_index()


df_squirls.to_parquet(
    snakemake.output['model_postprocess'], 
    index=False)