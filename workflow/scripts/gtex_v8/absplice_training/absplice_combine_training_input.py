import pandas as pd
from absplice.utils import get_abs_max_rows
import ast

index = ['variant', 'gene_id', 'sample','tissue']

# universe (gene, sample) with at least one rare variant
universe = pd.read_csv(snakemake.input['universe'])
universe['variants_on_gene'] = universe['variants_on_gene'].apply(lambda x: ast.literal_eval(x))
universe = universe.explode('variants_on_gene')
universe = universe.rename(columns={'variants_on_gene': 'variant'})
universe['tissue'] = snakemake.params['tissue']
universe = universe.set_index(index)

# outliers
outliers = pd.read_csv(snakemake.input['outliers'])
outliers['outlier'] = True
outliers['tissue'] = snakemake.wildcards['tissue']
outliers = outliers[['variant', 'gene_id', 'sample', 'tissue', 'outlier']]
outliers = outliers.drop_duplicates().set_index(index)

# predictions
try:
    absplice_input = pd.read_parquet(snakemake.input['absplice_input']).reset_index()
except:
    absplice_input = pd.read_csv(snakemake.input['absplice_input'])
absplice_input = absplice_input[
    absplice_input['tissue'] == snakemake.wildcards['tissue']]
absplice_input = absplice_input.set_index(index)

# left join everything to universe
df_absplice_input = universe.join(outliers, how='left')
df_absplice_input = df_absplice_input.join(absplice_input, how='left')
df_absplice_input = df_absplice_input.reset_index()

df_absplice_input.to_csv(snakemake.output['absplice_input'], index=False)