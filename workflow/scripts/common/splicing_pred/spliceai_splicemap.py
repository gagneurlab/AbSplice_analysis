from absplice_scripts.dataclasses.junction import get_splice_site_intervals, \
    get_unique_splice_site_intervals_in_event, intervals_to_pyranges, Junction
from splicemap import SpliceCountTable as CountTable
import pandas as pd
from typing import List, Union, Iterable, Iterator
import pyranges as pr
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors import variants_to_pyranges

from splicemap.splice_map import SpliceMap


# SpliceAI
df_spliceai = pd.read_parquet(snakemake.input['spliceai'])

if df_spliceai.shape[0] > 0:
    # annotate gene_id
    df_gene_map = pd.read_csv(snakemake.input['gene_map'], sep='\t')
    gene_mapping = dict(zip(df_gene_map['gene_name'], df_gene_map['gene_id']))
    df_spliceai['gene_id'] = df_spliceai['gene_name'].map(gene_mapping)
    
    # SpliceMap
    df_splicemap5 = SpliceMap.read_csv(snakemake.input['splicemap_5'][0]).df
    df_splicemap3 = SpliceMap.read_csv(snakemake.input['splicemap_3'][0]).df

    df_splicemap = pd.concat([
        df_splicemap5,
        df_splicemap3
    ])
    
    df_splicemap = df_splicemap[
        df_splicemap['gene_id'].isin(df_spliceai['gene_id'].unique())
    ]
    
    valid_junc_genes = set(df_splicemap.set_index(['junctions', 'gene_id']).index) # subset for those in end, to not have wrong genes on other strand
    
    df_splicemap['events'] = df_splicemap['events'].str.split(';')
    df_splicemap = df_splicemap.set_index('junctions')

    df_splicemap['interval'] = df_splicemap['events'].apply(lambda x: get_unique_splice_site_intervals_in_event(x, overhang=(100, 100)))
    df_intervals_all = df_splicemap.explode('interval')[['interval']]

    # get unique intervals into pyranges
    i_list = df_intervals_all['interval'].values
    i_list = list(set(i_list))
    pr_intervals_unique = intervals_to_pyranges(i_list)
    df_intervals_unique = pr_intervals_unique.df

    # get unique variants into pyranges
    v_list = df_spliceai['variant'].apply(lambda x: Variant.from_str(x)).values
    v_list = list(set(v_list))
    pr_vars = variants_to_pyranges(v_list)

    # get overlap of intervals and variants
    pr_intervals_with_vars = pr_intervals_unique.join(pr_vars)
    df_intervals_with_vars = pr_intervals_with_vars.df

else:
    df_intervals_with_vars = pd.DataFrame()


if df_intervals_with_vars.shape[0] > 0:
    # interval to string for joining
    df_intervals_all['interval'] = df_intervals_all['interval'].apply(lambda x: x.__str__())
    df_intervals_unique['interval'] = df_intervals_unique['interval'].apply(lambda x: x.__str__())
    df_intervals_with_vars['interval'] = df_intervals_with_vars['interval'].apply(lambda x: x.__str__())

    # get junction information from df_intervals_all
    df_intervals_all = df_intervals_all.reset_index().set_index('interval')
    df_intervals_with_vars = df_intervals_with_vars.set_index('interval')
    df_spliceai_splicemap = df_intervals_all.join(df_intervals_with_vars, how='inner')\
                                                .reset_index()[['junctions', 'variant']].drop_duplicates()

    # get variant informations from df_spliceai
    df_spliceai_splicemap['variant'] = df_spliceai_splicemap['variant'].apply(lambda x: x.__str__())
    df_spliceai_splicemap = df_spliceai_splicemap.set_index('variant').join(df_spliceai.set_index('variant')).reset_index()

    # combine splicemap5 and splicemap3
    df_splicemap5 = df_splicemap5.set_index(['junctions', 'gene_id'])[[
        'ref_psi', 'median_n',
    ]].add_prefix('psi5_')
    df_splicemap3 = df_splicemap3.set_index(['junctions', 'gene_id'])[[
        'ref_psi', 'median_n',
    ]].add_prefix('psi3_')
    df_splicemap_ref_psi = df_splicemap5.join(df_splicemap3)

    # join spliceai_splicemap with psi_ref info
    df_spliceai_splicemap = df_spliceai_splicemap.set_index(['junctions', 'gene_id'])
    df_spliceai_splicemap = df_spliceai_splicemap.join(df_splicemap_ref_psi)

    # annotate tissue
    df_spliceai_splicemap['tissue'] = snakemake.wildcards['tissue']

    # subset for valid junctions and genes
    valid_junc_genes = set(df_spliceai_splicemap.index.intersection(valid_junc_genes))
    df_spliceai_splicemap = df_spliceai_splicemap.loc[valid_junc_genes]
    df_spliceai_splicemap = df_spliceai_splicemap.reset_index()

else:
    df_spliceai_splicemap = pd.DataFrame()

df_spliceai_splicemap.to_parquet(snakemake.output['spliceai_splicemap'], index=False)
