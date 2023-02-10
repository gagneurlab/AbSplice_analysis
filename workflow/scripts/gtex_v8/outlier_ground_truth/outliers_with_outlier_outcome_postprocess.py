import pandas as pd
from tqdm import tqdm
tqdm.pandas()
import re
import numpy as np
import ast

def string_to_list(x):
    return [i.replace('[', '').replace(']', '').replace(' ', '').replace("'", '') for i in x.split(',')]

def find_match(df, categories):
    for c in categories:
        if c in df['aberrantSpliceType']:
            return True
    return False


outliers = pd.read_csv(snakemake.input['outliers'])
annotation = pd.read_csv(snakemake.input['outliers_with_annotation'])

# check that all outliers are there
index = ['variant', 'gene_id', 'sample']
index_outlier = set(outliers.set_index(index).index)
index_annotation = set(annotation.set_index(index).index)
assert len(index_outlier.difference(index_annotation)) == 0
assert len(index_annotation.difference(index_outlier)) >= 0

# preprocess annotation. string of list to list
annotation['aberrantSpliceType'] = annotation['aberrantSpliceType'].replace('[nan]', '["none"]')
annotation['aberrantSpliceType'] = annotation['aberrantSpliceType'].apply(lambda x: string_to_list(x))

# only use necessary columns
annotation = annotation[[
    'variant',
    'gene_id',
    'sample',
    'aberrantSpliceType',
    'abs_Distance'
]]

# annotate outlier outcome categories
outlier_categories = [
    # 'all',
    'singleExonSkipping',
    'exonSkipping', 
    'exonElongation', 
    'exonTruncation', 
    'theta_increase',
    'theta_decrease',
    'psi_increase',
    'psi_decrease'
]

if annotation.shape[0] > 0:
    for c in tqdm(outlier_categories):
        annotation.loc[:, f'{c}'] = annotation.apply(lambda df: find_match(df, [f'{c}']), axis=1)

    annotation['exonSkipping_all'] = False
    annotation.loc[
        ((annotation['singleExonSkipping'] == True) | (annotation['exonSkipping'] == True)),
        'exonSkipping_all'
    ] = True

    # save results
    assert annotation.shape[0] == len(set(annotation.set_index(index).index))
annotation.to_csv(snakemake.output['outliers_all'], index=False)