import pandas as pd

def get_type(df):
    change_types = ['acceptor_gain', 'acceptor_loss', 'donor_gain', 'donor_loss']
    change_type = list()
    for c in change_types:
        if df['delta_score'] == df[c]:
            change_type.append(c)
    return change_type  

def mask_ref_psi(df, ref_psi_tolerance=0.1):
    mask = False
    if len(df['change_type']) == 1:
        change_type = df['change_type'][0]
        if change_type in ['acceptor_gain', 'donor_gain']:
            if df['psi3_ref_psi'] >= (1-ref_psi_tolerance) and df['psi5_ref_psi'] >= (1-ref_psi_tolerance):
                mask = True
        elif change_type in ['acceptor_loss', 'donor_loss']:
            if df['psi3_ref_psi'] <= ref_psi_tolerance and df['psi5_ref_psi'] <= ref_psi_tolerance:
                mask = True
    return mask

# SpliceAI
df_spliceai_splicemap = pd.read_parquet(snakemake.input['spliceai_splicemap'])

# annotate change type of SpliceAI ('acceptor_gain', 'acceptor_loss', 'donor_gain', 'donor_loss')
df_spliceai_splicemap['change_type'] = df_spliceai_splicemap.apply(
    lambda x: get_type(x), 
    axis=1)

df_spliceai_splicemap['mask_ref_psi'] = df_spliceai_splicemap.apply(
    lambda df: mask_ref_psi(df, ref_psi_tolerance=float(snakemake.params['ref_psi_tolerance'])), 
    axis=1)

df_spliceai_splicemap_ref_psi = df_spliceai_splicemap[
    df_spliceai_splicemap['mask_ref_psi'] == False
]

df_spliceai_splicemap_ref_psi.to_parquet(snakemake.output['spliceai_splicemap_ref_psi'], index=False)