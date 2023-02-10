import pandas as pd
from tqdm import tqdm
tqdm.pandas()

def find_match(df, consequences):
    for c in consequences:
        if c in df['Consequence']:
            return True
    return False

# read universe
universe = pd.read_csv(snakemake.input['universe_variant'])\
    .rename(columns={'variants_on_gene': 'variant'})

# read VEP
df_vep = pd.read_parquet(snakemake.input['vep'])
df_vep['variant'] = snakemake.wildcards['vcf_id'] + ':' + df_vep['end'].astype(str) + ':' + df_vep['ref'] + '>' + df_vep['alt']
df_vep = df_vep[['variant', 'Consequence']]

# join universe with VEP consequence
universe = universe.set_index('variant').join(
    df_vep.set_index('variant')).reset_index()

# annotate variant categories
universe['Consequence'] = universe['Consequence'].fillna('missing')
cons_anno = [x for x in snakemake.params['var_categories'] if x != 'all']
for c in cons_anno:
    universe.loc[:, f'{c}'] = universe.progress_apply(lambda df: find_match(df, [f'{c}']), axis=1)

# filter based on variant category
consequences = snakemake.wildcards['var_category'].split('__')
print(consequences)
if snakemake.wildcards['var_category'] != 'all':
    print('filtering for consequences ...')
    universe = universe[
        universe.progress_apply(lambda df: find_match(df, consequences), axis=1)
    ]

# universe = universe.drop_duplicates()
# save universe category
universe.to_csv(snakemake.output['universe_vep'], index=False)
