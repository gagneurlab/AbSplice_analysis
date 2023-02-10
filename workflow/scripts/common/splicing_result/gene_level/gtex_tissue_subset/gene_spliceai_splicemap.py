import pandas as pd
from absplice import SplicingOutlierResult
from absplice.utils import get_abs_max_rows
from absplice_scripts.utils.mapping_utils import subset_tissues

df = pd.read_parquet(snakemake.input['pred_spliceai_splicemap'])
df_vcf_annotation = pd.read_csv(snakemake.input['var_samples_df'])

# join sample info
df = df.set_index('variant').join(
    df_vcf_annotation.set_index('variant'), 
    how='inner')\
    .reset_index()

tissue_map = pd.read_csv(snakemake.input['tissue_map'])
tissue_map = dict(zip(tissue_map['tissue'], tissue_map['tissue_main']))
df = subset_tissues(df, tissue_map, chosen_tissue=snakemake.wildcards['gtex_tissue'])

df = get_abs_max_rows(
    df = df.set_index(['gene_id', 'sample']),
    groupby = ['gene_id', 'sample'],
    max_col = 'delta_score',
)

assert df.shape[0] == len(df.index.unique())

df.to_parquet(snakemake.output['result_spliceai_splicemap_gene'])

