from absplice_scripts.dataclasses.junction import get_splice_site_intervals, \
    get_unique_splice_site_intervals_in_event, intervals_to_pyranges
from splicemap import SpliceCountTable as CountTable
import pandas as pd
from typing import List, Union, Iterable, Iterator
import pyranges as pr
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors import variants_to_pyranges


ct = CountTable.read_csv(snakemake.input['count_table_updated'])
df_rare_vars = pd.read_csv(snakemake.input['rare_vars'])
df_rare_vars = df_rare_vars[df_rare_vars['sample'].isin(ct.samples)]

rare_vars_chroms = set(df_rare_vars['variant'].apply(lambda x: x.split(':')[0]))
ct.df = ct.df[
    ct.df['Chromosome'].isin(rare_vars_chroms)]
    
df = pd.concat(
    [pd.DataFrame(ct.event5['events'].str.split(';')), 
     pd.DataFrame(ct.event3['events'].str.split(';'))]
)

df['interval'] = df['events'].apply(lambda x: get_unique_splice_site_intervals_in_event(x, overhang=(100, 100)))
df_intervals_all = df.explode('interval')[['interval']]

# get unique intervals into pyranges
i_list = df_intervals_all['interval'].values
i_list = list(set(i_list))
pr_intervals_unique = intervals_to_pyranges(i_list)
df_intervals_unique = pr_intervals_unique.df

# get unique variants into pyranges
v_list = df_rare_vars['variant'].apply(lambda x: Variant.from_str(x)).values
v_list = list(set(v_list))
pr_rare_vars = variants_to_pyranges(v_list)

# get overlap of intervals and variants
pr_intervals_with_rare = pr_intervals_unique.join(pr_rare_vars)
df_intervals_with_rare = pr_intervals_with_rare.df

# interval to string for joining
df_intervals_all['interval'] = df_intervals_all['interval'].apply(lambda x: x.__str__())
df_intervals_unique['interval'] = df_intervals_unique['interval'].apply(lambda x: x.__str__())
df_intervals_with_rare['interval'] = df_intervals_with_rare['interval'].apply(lambda x: x.__str__())

# get junction information from df_intervals_all
df_intervals_all = df_intervals_all.reset_index().set_index('interval')
df_intervals_with_rare = df_intervals_with_rare.set_index('interval')
df_junctions_with_rare_event = df_intervals_all.join(df_intervals_with_rare, how='inner')\
                                            .reset_index()[['junctions', 'variant']].drop_duplicates()

# get variant informations from df_rare_vars
df_junctions_with_rare_event['variant'] = df_junctions_with_rare_event['variant'].apply(lambda x: x.__str__())
df_junctions_with_rare_event = df_junctions_with_rare_event.set_index('variant').join(df_rare_vars.set_index('variant')).reset_index()
df_junctions_with_rare_event.to_csv(snakemake.output['junctions_with_rare_variant'], index=False)
