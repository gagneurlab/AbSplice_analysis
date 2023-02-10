from tqdm import tqdm
import pandas as pd
from splicemap import SpliceCountTable
from absplice_scripts.utils.variant_utils import filter_long_variants
from mmsplice.utils import left_normalized
from kipoiseq import Variant

df_rare_vars = pd.read_csv(snakemake.input['rare_variants'])
ct = SpliceCountTable.read_csv(snakemake.input['count_table'])

# filter for samples that have DNA and RNA
df_rare_vars['sample'] = df_rare_vars['sample'].astype(str)
samples_DNA = set(df_rare_vars['sample'])
samples_RNA = set(ct.samples)
samples = samples_DNA.intersection(samples_RNA)
df_rare_vars = df_rare_vars[df_rare_vars['sample'].isin(samples)]

df_rare_vars.to_csv(snakemake.output['rare_vars_tissue'], index=False)

# construct benchmark universe
df_universe = df_rare_vars[['variant', 'gene_id', 'sample']]
df_universe = df_universe.groupby(['gene_id', 'sample'])['variant'].apply(lambda x: sorted(set(x)))
df_universe = pd.DataFrame(df_universe)
df_universe = df_universe.rename(columns={'variant': 'variants_on_gene'})
df_universe.to_csv(snakemake.output['universe'])