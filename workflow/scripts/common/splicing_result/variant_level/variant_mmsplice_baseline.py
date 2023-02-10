import pandas as pd
from tqdm import tqdm
import pandas as pd
from absplice.utils import get_abs_max_rows

df_pred = pd.read_csv(snakemake.input['pred_mmsplice_baseline'])\
            .rename(columns={'ID': 'variant'})\
            .set_index('variant')
df_vcf_annotation = pd.read_csv(snakemake.input['var_samples_df'])\
            .set_index('variant')

df_pred = df_pred.join(df_vcf_annotation, how='inner').reset_index()

# clean up gene_ids
df_pred['gene_id'] = df_pred['gene_id'].apply(lambda x: x.split(';'))
df_pred = df_pred.explode('gene_id')
df_pred['gene_id'] = df_pred['gene_id'].apply(lambda x: x.split('.')[0])

df_pred = get_abs_max_rows(
    df = df_pred.set_index(['variant', 'gene_id', 'sample']),
    groupby = ['variant', 'gene_id', 'sample'],
    max_col = 'delta_logit_psi',
#     dropna=False
)

assert df_pred.shape[0] == len(df_pred.index.unique())

df_pred.to_csv(snakemake.output['result_mmsplice_baseline_variant'])
