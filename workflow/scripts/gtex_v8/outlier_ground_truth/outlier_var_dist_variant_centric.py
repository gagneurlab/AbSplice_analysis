import pandas as pd
from tqdm import tqdm
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors import variants_to_pyranges
from absplice_scripts.dataclasses.junction import get_splice_site_intervals, \
    get_unique_splice_site_intervals_in_event, intervals_to_pyranges

df_outlier = pd.read_csv(snakemake.input['outliers_signif'])
df_outlier['junctions'] = df_outlier['junctions_j']
df_outlier['events'] = df_outlier['events_j']
df_outlier['events'] = df_outlier['events'].str.split(';')
df_outlier = df_outlier[
    df_outlier['seqnames_j'] == snakemake.wildcards['vcf_id']
]
df_outlier = df_outlier.set_index('junctions')

df_rare_vars = pd.read_csv(snakemake.input['rare_vars'])

common_samples = sorted(set(df_outlier['sample']).intersection(
    set(df_rare_vars['sample'])))


def get_abs_min_rows(df, groupby, min_col, dropna=True):
    df = df.reset_index()
    _df = df.copy()
    _df[min_col] = _df[min_col].abs()
    min_scores = _df.groupby(groupby, dropna=dropna)[min_col].idxmin()
    return df.iloc[min_scores.values].set_index(groupby)


def junc_var_dist(df_outlier, df_rare_vars, sample, k=100):
    df_outlier = df_outlier[df_outlier['sample'] == sample]
    df_rare_vars = df_rare_vars[df_rare_vars['sample'] == sample]
    
    # get unique variants into pyranges
    v_list = df_rare_vars['variant'].apply(lambda x: Variant.from_str(x)).values
    v_list = list(set(v_list))
    pr_rare_vars = variants_to_pyranges(v_list)
    
    # get unique intervals of outlier events
    df_outlier['interval'] = df_outlier['events'].apply(lambda x: get_unique_splice_site_intervals_in_event(x, overhang=(0, 0)))
    df_intervals_all = df_outlier.explode('interval')[['interval']]
    
    # get unique intervals into pyranges
    i_list = df_intervals_all['interval'].values
    i_list = list(set(i_list))
    pr_intervals_unique = intervals_to_pyranges(i_list)
    df_intervals_unique = pr_intervals_unique.df
    
    # get nearest distance of junctions to rare variants
    pr_intervals_with_rare = pr_rare_vars.k_nearest(pr_intervals_unique, k=k)
    df_intervals_with_rare = pr_intervals_with_rare.df
    
    # interval to string for joining
    df_intervals_with_rare['variant'] = df_intervals_with_rare['variant'].apply(lambda x: x.__str__())
    df_intervals_all['interval'] = df_intervals_all['interval'].apply(lambda x: x.__str__())
    df_intervals_unique['interval'] = df_intervals_unique['interval'].apply(lambda x: x.__str__())
    df_intervals_with_rare['interval'] = df_intervals_with_rare['interval'].apply(lambda x: x.__str__())
    
    # get junction information from df_intervals_all
    df_intervals_all = df_intervals_all.reset_index().set_index('interval')
    df_intervals_with_rare = df_intervals_with_rare.set_index('interval')
    df_junctions_with_rare_event = df_intervals_all.join(
        df_intervals_with_rare, how='inner').reset_index()
    
    # get variant informations from df_rare_vars
    df_junctions_with_rare_event['variant'] = df_junctions_with_rare_event['variant'].apply(lambda x: x.__str__())
    df_junctions_with_rare_event = df_junctions_with_rare_event.set_index('variant').join(
        df_rare_vars.set_index('variant')).reset_index()
    
    # get minimum distance of variant and junction
    df_junctions_with_rare_event = get_abs_min_rows(
        df_junctions_with_rare_event.set_index(['sample', 'variant', 'junctions']),
        min_col='Distance',
        groupby=['sample', 'variant', 'junctions']
    ).reset_index()
        
    return df_junctions_with_rare_event


success = False
failure = False
sample_success = False
k_range = [500, 1000, 2000, 3000, 10000]
count = 0

while failure == False and success == False:
    
    df_junc_vars = list()
    for sample in tqdm(common_samples):
        
        sample_success = False
        count = 0
        
        while failure == False and sample_success == False:
            k = k_range[count]
            res_temp = junc_var_dist(df_outlier, df_rare_vars, sample, k)
        
            try:
                assert len(set(df_outlier[df_outlier['sample'] == sample].reset_index()['junctions']).difference(
                    set(res_temp['junctions']))) == 0
                assert len(set(df_rare_vars[df_rare_vars['sample'] == sample].reset_index()['variant']).difference(
                    set(res_temp['variant']))) == 0
                df_junc_vars.append(res_temp)
                sample_success = True
                print(f'success: {sample}, k={k}')
            except:
                print(f'failed: {sample}, k={k}')
                count += 1
                if count < len(k_range):
                    continue
                else:
                    failure = True
                    break
    df_junc_vars = pd.concat(df_junc_vars)
    
    if len(set(common_samples).difference(set(df_junc_vars['sample']))) == 0:
        success = True
    print(f'success: {success}')

assert success == True


def get_abs_min_rows(df, groupby, min_col, dropna=True):
    df = df.reset_index()
    _df = df.copy()
    _df[min_col] = _df[min_col].abs()
    min_scores = _df.groupby(groupby, dropna=dropna)[min_col].idxmin()
    return df.iloc[min_scores.values].set_index(groupby)


df_junc_vars_min = get_abs_min_rows(
    df_junc_vars.set_index(['sample', 'variant', 'junctions']),
    min_col='Distance',
    groupby=['sample', 'variant', 'junctions']
).reset_index()


df_junc_vars_min['abs_Distance'] = abs(df_junc_vars_min['Distance'])

df_junc_vars_min = df_junc_vars_min.set_index(['junctions', 'sample']).join(
    df_outlier.reset_index()[['junctions', 'sample', 'gene_id']].drop_duplicates().set_index(['junctions', 'sample'])[['gene_id']]).reset_index()


df_junc_vars_min.to_csv(snakemake.output['variant_outlier_dist'], index=False)