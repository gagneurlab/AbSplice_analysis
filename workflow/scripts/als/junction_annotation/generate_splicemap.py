import pandas as pd
import pyranges as pr
from splicemap import SpliceCountTable as CountTable

if snakemake.params['method'] == 'kn':
    used_method = 'k/n'
else:
    used_method = snakemake.params['method']

df_gene_expression = pd.read_csv(snakemake.input['gene_expression'])
df_gene_expression = df_gene_expression[['gene_id', snakemake.params['tissue']]]

ct = CountTable.read_csv(snakemake.input['count_table'], name=snakemake.params['tissue'])
ct = CountTable(
    ct.df[['Chromosome', 'Start', 'End', 'Strand']
          + [i for i in ct.samples if i.startswith('CTRL')]], 
    name=snakemake.params['tissue'],
    gene_expression=df_gene_expression)

if snakemake.params['event_filter'] == 'median_cutoff':
    ct_psi5 = ct.event5_median_filter(
        cutoff=snakemake.params['median_cutoff'])
    ct_psi3 = ct.event3_median_filter(
        cutoff=snakemake.params['median_cutoff'])
elif snakemake.params['event_filter'] == 'gaussian_mixture':
    [ct_psi5, cutoff_5] = ct.event5_count_filter()
    [ct_psi3, cutoff_3] = ct.event3_count_filter()

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
