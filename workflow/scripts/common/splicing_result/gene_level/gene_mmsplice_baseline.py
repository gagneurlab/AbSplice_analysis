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

# Aggregate predictions
df_pred_agg = get_abs_max_rows(df_pred.set_index(['gene_id', 'sample']), 
                               groupby=['gene_id', 'sample'], 
                               max_col='delta_logit_psi')

df_pred_agg.to_parquet(snakemake.output['result_mmsplice_baseline_gene'])
