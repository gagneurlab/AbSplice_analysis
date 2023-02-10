import pandas as pd
from tqdm import tqdm
import pandas as pd
from absplice.utils import get_abs_max_rows

df_pred = pd.read_parquet(snakemake.input['pred_spliceai_splicemap'])
df_vcf_annotation = pd.read_csv(snakemake.input['var_samples_df'])

# join sample info
df_pred = df_pred.set_index('variant').join(
    df_vcf_annotation.set_index('variant'), 
    how='inner')\
    .reset_index()

# Aggregate predictions
df_pred_agg = get_abs_max_rows(
    df_pred.set_index(['variant', 'gene_id', 'sample', 'tissue']), 
    groupby=['variant', 'gene_id', 'sample', 'tissue'], 
    max_col='delta_score'
)

df_pred_agg.to_csv(snakemake.output['result_spliceai_splicemap_variant'])
