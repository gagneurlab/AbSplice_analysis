from kipoiseq.extractors.vcf_seq import MultiSampleVCF
from gnomad_rocksdb import GnomadMafDB
from absplice_scripts.data.vcf_utils import get_variant_around_event_psi5, \
    get_variant_around_event_psi3, VariantMafFilter, PrivateVariantFilter, ReadDepthFilter, GQFilter

maf = GnomadMafDB(snakemake.input['maf'])
vcf = MultiSampleVCF(snakemake.input['vcf'])

vcf.query_all() \
   .filter(lambda v: v.filter is None) \
   .filter(lambda v: VariantMafFilter(cutoff=0.001, population=maf)(v)) \
   .filter(lambda v: PrivateVariantFilter(vcf, max_num_samples=2)(v)) \
   .to_sample_csv(snakemake.output['vcf_annotation'], format_fields=snakemake.params['format_fields'])
