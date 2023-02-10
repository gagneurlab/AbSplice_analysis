import pandas as pd
import pyranges as pr
import re

def extract_dna_id(filename):
    pattern = "genomics\/3_vcf\/.*\/(.*)\/.*"
    dna_id = re.search(pattern, filename).group(1)
    return dna_id

df_anno = pd.read_csv(snakemake.input['sample_anno'], sep='\t')

_df_anno = df_anno[
    ~df_anno['DNA_VCF_FILE'].isna()
]
_df_anno['VCF_ID'] = _df_anno.apply(lambda x: extract_dna_id(x['DNA_VCF_FILE']), axis=1)
df_anno = df_anno.set_index('RNA_ID').join(_df_anno.set_index('RNA_ID')['VCF_ID']).reset_index()
df_anno.to_csv(snakemake.output['sample_anno_updated'], sep='\t')
