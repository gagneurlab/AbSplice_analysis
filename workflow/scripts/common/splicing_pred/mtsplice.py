# import pyarrow
import kipoi.metadata
from mmsplice import predict_save, MMSplice
from mmsplice.vcf_dataloader import SplicingVCFDataloader

dl = SplicingVCFDataloader(snakemake.input['gtf'], snakemake.input['fasta'],
                           snakemake.input['vcf'], tissue_specific=True)
model = MMSplice()
predict_save(model, dl, snakemake.output['result'])
