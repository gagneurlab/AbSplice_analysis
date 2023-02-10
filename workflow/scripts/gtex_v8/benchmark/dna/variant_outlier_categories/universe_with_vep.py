import pandas as pd
from tqdm import tqdm
tqdm.pandas()

# VEP consequences ranked for most severe
consequence_ranking = {
    'transcript_ablation': 1,
    'splice_acceptor_variant': 2,
    'splice_donor_variant': 3,
    'stop_gained': 4,
    'frameshift_variant': 5,
    'stop_lost': 6,
    'start_lost': 7,
    'transcript_amplification': 8,
    'inframe_insertion': 9,
    'inframe_deletion': 10,
    'missense_variant': 11,
    'protein_altering_variant': 12,
    'splice_region_variant': 13,
    'splice_donor_5th_base_variant': 14,
    'splice_donor_region_variant': 15,
    'splice_polypyrimidine_tract_variant': 16,
    'incomplete_terminal_codon_variant': 17,
    'start_retained_variant': 18,
    'stop_retained_variant': 19,
    'synonymous_variant': 20,
    'coding_sequence_variant': 21,
    'mature_miRNA_variant': 22,
    '5_prime_UTR_variant': 23,
    '3_prime_UTR_variant': 24,
    'non_coding_transcript_exon_variant': 25,
    'intron_variant': 26,
    'NMD_transcript_variant': 27,
    'non_coding_transcript_variant': 28,
    'upstream_gene_variant': 29,
    'downstream_gene_variant': 30,
    'TFBS_ablation': 31,
    'TFBS_amplification': 32,
    'TF_binding_site_variant': 33,
    'regulatory_region_ablation': 34,
    'regulatory_region_amplification': 35,
    'feature_elongation': 36,
    'regulatory_region_variant': 37,
    'feature_truncation': 38,
    'intergenic_varian': 39
}
consequence_ranking_inv = {i:c for c,i in consequence_ranking.items()}

# find the most severe consequence
def find_most_severe(df):
    consequences = df['Consequence']
    min_rank = 10000
    for c in consequences:
        temp_rank = consequence_ranking[c]
        if temp_rank < min_rank:
            min_rank = temp_rank
    min_consequence = consequence_ranking_inv[min_rank]
    return min_consequence

# from all consequences check if given consequence included
def find_match(df, consequences):
    for c in consequences:
        if c in df['Consequence']:
            return True
    return False

# from all consequences check if given consequence is only one annotated
def find_pure_match(df, consequence):
    if len(df['Consequence']) == 1:
        if df['Consequence'][0] == consequence:
            return True
    return False

def is_consequence_most_severe(df, consequence):
    # __import__("pdb").set_trace()
    return consequence == df['most_severe_consequence']


# read universe
universe = pd.read_csv(snakemake.input['universe_variant'])

universe_dim = universe.shape[0]

# read VEP
df_vep = pd.read_parquet(snakemake.input['vep'])
df_vep['variant'] = snakemake.wildcards['vcf_id'] + ':' + df_vep['end'].astype(str) + ':' + df_vep['ref'] + '>' + df_vep['alt']
#subset for protein coding
df_vep = df_vep[df_vep['BIOTYPE'] == 'protein_coding']
#get unique consequences for each variant
df_vep = df_vep[['variant', 'Consequence']].explode('Consequence')
df_vep = df_vep.groupby('variant').apply(lambda df: sorted(set(df['Consequence'])))
df_vep = pd.DataFrame(df_vep).rename(columns={0: 'Consequence'}).reset_index()

df_vep = df_vep[['variant', 'Consequence']]

df_vep['most_severe_consequence'] = df_vep.apply(lambda df: find_most_severe(df), axis=1)

# annotate variant category
cons_anno = [
    'intron_variant', 
    'splice_region_variant',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'stop_lost',
    'missense_variant',
    'synonymous_variant',
]
for c in tqdm(cons_anno):
    df_vep.loc[:, f'contains_{c}'] = df_vep.apply(lambda df: find_match(df, [f'{c}']), axis=1)
    df_vep.loc[:, f'pure_{c}'] = df_vep.apply(lambda df: find_pure_match(df, [f'{c}']), axis=1)
    df_vep.loc[:, f'{c}'] = df_vep.apply(lambda df: is_consequence_most_severe(df, f'{c}'), axis=1)
    
df_vep['exonic_variant'] = False
df_vep.loc[
    ((df_vep['stop_gained'] == True) | (df_vep['stop_lost'] == True) | (df_vep['missense_variant'] == True) | (df_vep['synonymous_variant'] == True)),
    'exonic_variant'
] = True

# join universe with VEP consequence
universe = universe.set_index('variant')
df_vep = df_vep.set_index('variant')

assert df_vep.shape[0] == len(set(df_vep.index))

universe = universe.join(df_vep).reset_index()

# save universe category
universe = universe.drop(columns='Consequence')
assert universe.shape[0] == universe_dim
universe.to_csv(snakemake.output['universe_vep'], index=False)