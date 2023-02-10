from tqdm import tqdm
import pandas as pd
import ast

universe = pd.read_csv(snakemake.input['universe'])
universe['variants_on_gene'] = universe['variants_on_gene'].apply(lambda x: ast.literal_eval(x))
universe = universe.explode('variants_on_gene')
universe = universe.rename(columns={'variants_on_gene': 'variant'})
universe.to_csv(snakemake.output['universe'], index=False)
