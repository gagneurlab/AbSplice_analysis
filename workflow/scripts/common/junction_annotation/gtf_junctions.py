from splicemap.count_table import SpliceCountTable as CountTable
import pyranges as pr
import pandas as pd

df = pd.DataFrame({
        'Chromosome': ['chr1', 'chr1'],
        'Start': [5, 22],
        'End': [30, 30],
        'Strand': ['+', '+'],
        's1': [1, 1],
        's2': [2, 1],
    })

ct = CountTable(df, name='dummy_count_table')

gr_gtf = pr.read_gtf(snakemake.input['gtf'])
df_gtf_junc = ct._load_junction_from_gtf(gr_gtf)

df_gtf_junc.to_csv(snakemake.output['gtf_junctions'])
