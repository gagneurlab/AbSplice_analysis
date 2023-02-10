import pandas as pd
from splicing_outlier_prediction.utils import get_abs_max_rows
import ast

df_absplice_input = pd.read_csv(snakemake.input['absplice_input'])

# left join rest of the models
squirls = pd.read_csv(snakemake.input['squirls'])
cadd_splice = pd.read_csv(snakemake.input['cadd_splice'])
mtsplice = pd.read_csv(snakemake.input['mtsplice'])
mtsplice = mtsplice[
    mtsplice['tissue'] == snakemake.wildcards['tissue']
]
spliceai_splicemap = pd.read_csv(snakemake.input['spliceai_splicemap'])
spliceai_splicemap = spliceai_splicemap[
    spliceai_splicemap['tissue'] == snakemake.wildcards['tissue']
]
spliceai_splicemap = spliceai_splicemap.rename(columns={'delta_score': 'delta_score_splicemap'})

spliceai_splicemap_ref_psi = pd.read_csv(snakemake.input['spliceai_splicemap_ref_psi'])
spliceai_splicemap_ref_psi = spliceai_splicemap_ref_psi[
    spliceai_splicemap_ref_psi['tissue'] == snakemake.wildcards['tissue']
]
spliceai_splicemap_ref_psi = spliceai_splicemap_ref_psi.rename(columns={'delta_score': 'delta_score_splicemap_ref_psi'})

index = ['variant', 'gene_id', 'sample']

df_absplice_input = df_absplice_input.set_index(index).join(
    squirls.set_index(index)[['squirls_scores']], how='left').join(
    mtsplice.set_index(index)[['delta_logit_psi_mtsplice']], how='left').join(
    cadd_splice.set_index(index)[['PHRED']], how='left').join(
    spliceai_splicemap.set_index(index)[['delta_score_splicemap']], how='left').join(
    spliceai_splicemap_ref_psi.set_index(index)[['delta_score_splicemap_ref_psi']], how='left').reset_index()

df_absplice_input.to_csv(snakemake.output['absplice_input'], index=False)
