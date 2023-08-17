import pyarrow
from mmsplice import predict_save, MMSplice
from mmsplice.vcf_dataloader import SplicingVCFDataloader
import pandas as pd

dl = SplicingVCFDataloader(snakemake.input['gtf'], snakemake.input['fasta'],
                           snakemake.input['vcf'])
model = MMSplice()
predict_save(model, dl, output_path=snakemake.output['result'])
