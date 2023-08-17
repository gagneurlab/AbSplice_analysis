from splicemap import SpliceCountTable as CountTable
from absplice_scripts.data.DROP_annotations import sample_individual_mapping

# Read and convert to pyranges
ct = CountTable.read_csv(snakemake.input['raw_count_table'],
                         name=snakemake.params['tissue'])

if snakemake.params['update_samples']:
    ct.update_samples(sample_individual_mapping(
        annotation_table_path=snakemake.params['annotation_table'],
        key_assay=snakemake.params['key_assay'],                           
        value_assay=snakemake.params['value_assay']))

# if snakemake.params['remove_chr_from_chrom_annotation']:
#     ct.df['Chromosome'] = ct.df['Chromosome'].str.replace('chr', '')
    
if 'chr' not in ct.df['Chromosome'].values[0]:
    ct.df['Chromosome'] = 'chr' + ct.df['Chromosome'].astype(str)

if snakemake.params['infer_strand']:
    ct.infer_strand(snakemake.input['fasta'])

ct.to_csv(snakemake.output['updated_count_table'])
