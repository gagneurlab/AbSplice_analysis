import pdb
import pandas as pd
import pyarrow
from absplice.result import SplicingOutlierResult
from absplice.cat_dataloader import CatInference

# Infer CAT
cat_dl = CatInference(
    splicemap5=list(snakemake.input['splicemap_5']),
    splicemap3=list(snakemake.input['splicemap_3']),
    count_cat=snakemake.input['cat_count_table'],
    name=snakemake.params['tissue_cat']
)

result = SplicingOutlierResult(
    df_mmsplice = snakemake.input['mmsplice_splicemap'],
)
result.add_samples(pd.read_csv(snakemake.input['var_samples_df']))

result.infer_cat(cat_dl, progress=True)
result.df_mmsplice_cat.to_parquet(snakemake.output['result'], partition_cols='tissue_cat', engine='pyarrow')

