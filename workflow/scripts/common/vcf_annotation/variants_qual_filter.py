from tqdm import tqdm
import pandas as pd
from absplice_scripts.utils.variant_utils import filter_long_variants
from mmsplice.utils import left_normalized
from kipoiseq import Variant

df_rare_vars = pd.read_csv(snakemake.input['rare_variants'])

# quality filter variants
if 'GQ' in df_rare_vars.columns:
    df_rare_vars = df_rare_vars[df_rare_vars['GQ'] >= snakemake.params['min_GQ']]
if 'AD' in df_rare_vars.columns:
    df_rare_vars['DP_ALT'] = df_rare_vars['AD'].apply(lambda x: x.split(',')[1]).astype(int)
    df_rare_vars = df_rare_vars[df_rare_vars['DP_ALT'] >= snakemake.params['min_DP_ALT']]

# filter long variants
df_rare_vars['left_normalized_variant'] = df_rare_vars['variant'].apply(
    lambda x: left_normalized(Variant.from_str(x)).__str__())
df_rare_vars = df_rare_vars[df_rare_vars['left_normalized_variant'].map(
    lambda v: filter_long_variants(v, max_length=snakemake.params['max_length_indel']))]

df_rare_vars.to_csv(snakemake.output['rare_variants_filtered'], index=False)