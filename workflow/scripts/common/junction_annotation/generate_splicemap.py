import pandas as pd
from splicemap import SpliceCountTable as CountTable

if snakemake.params['method'] == 'kn':
    used_method = 'k/n'
else:
    used_method = snakemake.params['method']
    
df_gene_expression = pd.read_csv(snakemake.input['gene_expression'])
df_gene_expression = df_gene_expression[['gene_id', snakemake.params['tissue']]]

ct_test = CountTable.read_csv(snakemake.input['count_table_updated'], 
                              name=snakemake.params['tissue'])
ct_test = CountTable(ct_test.df, ct_test.name, df_gene_expression)

df_test = pd.read_csv(snakemake.input['train_test_split'], index_col=0)
df_counts = ct_test.counts
df_counts[df_test == 1] = 0
ct_test.df[ct_test.samples] = df_counts

if snakemake.params['event_filter'] == 'median_cutoff':
    ct_psi5 = ct_test.event5_median_filter(
        cutoff=snakemake.params['median_cutoff'])
    ct_psi3 = ct_test.event3_median_filter(
        cutoff=snakemake.params['median_cutoff'])
elif snakemake.params['event_filter'] == 'gaussian_mixture':
    [ct_psi5, cutoff_5] = ct_test.event5_count_filter()
    [ct_psi3, cutoff_3] = ct_test.event3_count_filter()

ct_psi5 = ct_psi5.quantile_filter(quantile=snakemake.params['percentile'],
                                  min_read=snakemake.params['percentile_min_read'])
ct_psi3 = ct_psi3.quantile_filter(quantile=snakemake.params['percentile'],
                                  min_read=snakemake.params['percentile_min_read'])

gtf_file = snakemake.input['gtf_file']
ct_psi5.infer_annotation(gtf_file)
ct_psi3.infer_annotation(gtf_file)

df_psi5 = ct_psi5.ref_psi5(method=used_method)
df_psi3 = ct_psi3.ref_psi3(method=used_method)

df_psi5.to_csv(snakemake.output['splicemap_psi5'])
df_psi3.to_csv(snakemake.output['splicemap_psi3'])
