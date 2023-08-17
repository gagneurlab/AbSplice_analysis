from kipoiseq import Variant
import pandas as pd

def filter_long_variants(variant, max_length=10):
    if pd.isna(variant):
        return True
    
    variant = Variant.from_str(variant)
    length = max(len(variant.ref), len(variant.alt))
    return length < max_length


chromlenghts = {}
chromlenghts['hg19'] = {}
chromlenghts['hg38'] = {}
chromlenghts['hg38'] = {
    '1': 248956422,
    '2': 242193529,
    '3': 198295559,
    '4': 190214555,
    '5': 181538259,
    '6': 170805979,
    '7': 159345973,
    '8': 145138636,
    '9': 138394717,
    '10': 133797422,
    '11': 135086622,
    '12': 133275309,
    '13': 114364328,
    '14': 107043718,
    '15': 101991189,
    '16': 90338345,
    '17': 83257441,
    '18': 80373285,
    '19': 58617616,
    '20': 64444167,
    '21': 46709983,
    '22': 50818468,
    'X': 156040895,
    'Y': 57227415,
}

chromlenghts['hg19'] = {
    '1': 249250621,
    '2': 243199373,
    '3': 198022430,
    '4': 191154276,
    '5': 180915260,
    '6': 171115067,
    '7': 159138663,
    '8': 146364022,
    '9': 141213431,
    '10': 135534747,
    '11': 135006516,
    '12': 133851895,
    '13': 115169878,
    '14': 107349540,
    '15': 102531392,
    '16': 90354753,
    '17': 81195210,
    '18': 78077248,
    '19': 59128983,
    '20': 63025520,
    '21': 48129895,
    '22': 51304566,
    'X': 155270560,
    'Y': 59373566,
}

def var(row, chrom_col='chrom', pos_col='pos', ref_col='ref', alt_col='alt'):
    return row[chrom_col] + ':' + str(row[pos_col]) + ':' + str(row[ref_col]) + '>' + str(row[alt_col])

def var_from_str(row):
    return Variant.from_str(row['variant'])

def to_vcf(variants, filepath, chr_annotation, genome_version):
    
    with open(filepath, 'w') as vcf_file:
        vcf_file.write('##fileformat=VCFv4.0\n')
        for chrom, chromlength in chromlenghts[genome_version].items():
            vcf_file.write(f'##contig=<ID={chr_annotation}{chrom},length={chromlength}>\n')
        vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        for v in variants:
            # v = chr_annotation + v
            vcf_file.write('\t'.join([
                v.chrom, str(v.pos), '.', v.ref, v.alt, '.', '.', '.'
            ]))
            vcf_file.write('\n')