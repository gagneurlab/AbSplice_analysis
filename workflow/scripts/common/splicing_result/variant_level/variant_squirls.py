import pandas as pd
from tqdm import tqdm
import pandas as pd
from absplice.utils import get_abs_max_rows

df_pred = pd.read_csv(snakemake.input['pred_squirls'])\
            .set_index('variant')
df_vcf_annotation = pd.read_csv(snakemake.input['var_samples_df'])\
            .set_index('variant')

df_pred = df_pred.join(df_vcf_annotation, how='inner').reset_index()

# Aggregate predictions
df_pred_agg = get_abs_max_rows(df_pred.set_index(['variant', 'gene_id', 'sample']), 
                               groupby=['variant', 'gene_id', 'sample'], 
                               max_col='squirls_scores')

df_pred_agg.to_csv(snakemake.output['result_squirls_variant'])
