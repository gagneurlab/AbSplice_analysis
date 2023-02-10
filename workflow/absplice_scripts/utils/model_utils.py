model_dict_dna = {
    'CADD-Splice': 'PHRED_cadd_splice',
    'SQUIRLS': 'squirls_scores_squirls',
    'MTSplice': 'delta_logit_psi_mtsplice_mtsplice',

    'SpliceAI': 'delta_score_spliceai',
    'SpliceAI SpliceMap': 'delta_score_spliceai_splicemap',
    'SpliceAI SpliceMap + Ref PSI': 'delta_score_spliceai_splicemap_ref_psi',
    
    'MMSplice': 'delta_logit_psi_mmsplice',
    'MMSplice SpliceMap': 'delta_logit_psi_mmsplice_splicemap',
    'MMSplice SpliceMap + Ref_PSI': 'delta_psi_mmsplice_splicemap_ref_psi',
    
    'AbSplice-DNA': 'AbSplice_DNA_absplice_dna',
}

models_rna_only = {
    'CAT pval': 'pValueGene_g_minus_log10_cat_pval',
    'CAT infer SpliceMap': 'delta_psi_cat_mmsplice_splicemap_cat',
    'AbSplice-RNA (single cat)': 'AbSplice_RNA_absplice_rna_single_cat',
    'AbSplice-RNA (all cats)': 'AbSplice_RNA_absplice_rna_all_cats',
}

model_dict_rna = {**model_dict_dna,**models_rna_only}

pred_tools_dna = [
    'spliceai',
    'spliceai_splicemap',
    'spliceai_splicemap_ref_psi',
    'cadd_splice',
    'squirls',
    'mtsplice',
    'mmsplice',
    'mmsplice_splicemap',
    'mmsplice_splicemap_ref_psi',
    'absplice_dna', 
]
pred_tools_rna_only = [
    'cat_pval',
    'mmsplice_splicemap_cat',
    'absplice_rna_single_cat',
    'absplice_rna_all_cats',
]
pred_tools_rna = [
    *pred_tools_dna,
    *pred_tools_rna_only
]

tissue_specific_tools = [
    'mtsplice',
    'mmsplice_splicemap',
    'mmsplice_splicemap_ref_psi',
    'spliceai_splicemap',
    'spliceai_splicemap_ref_psi',
    'absplice_dna'
]

replace_dict_variants = {
    'exonic_variant': 'exonic',
    'splice_donor_variant': 'splice donor',
    'splice_acceptor_variant': 'splice acceptor',
    'splice_region_variant': 'splice region',
    'intron_variant': 'intron',
    'pure_intron_variant': 'pure intron',
    'stop_gained': 'stop gained',
    'stop_lost': 'stop lost',
    'synonymous_variant': 'synonymous',
    'all': 'all',
}

replace_dict_outliers = {
    'exonElongation': 'Exon elongation',
    'exonSkipping_all': 'Exon skipping',
    'exonTruncation': 'Exon truncation',
    'theta_increase': 'Splicing efficiency increase',
    'theta_decrease': 'Splicing efficiency decrease',
    'all': 'all',
    'splicing_efficiency_change': 'splicing_efficiency_change',
    'any_psi5_psi3_change': 'any_psi5_psi3_change',
}


model_cutoffs = {
    'delta_score_spliceai': {
        'high': 0.8,
        'medium': 0.5,
        'low': 0.2,
    },
    'delta_logit_psi_mmsplice': {
        'high': 2,
        'medium': 1.5,
        'low': 1,
    },
    'PHRED_cadd_splice': {
        'high': 50,
        'medium': 20,
        'low': 10,
    }
}
model_cutoffs['delta_psi_mmsplice_splicemap_ref_psi'] = {
        'high': 0.2,
        'medium': 0.1,
        'low': 0.05,
}
model_cutoffs['AbSplice_DNA_absplice_dna'] = {
        'high': 0.2,
        'medium': 0.05,
        'low': 0.01,
}
model_cutoffs['delta_score_spliceai_splicemap_ref_psi'] = model_cutoffs['delta_score_spliceai']
model_cutoffs['delta_score_spliceai_splicemap'] = model_cutoffs['delta_score_spliceai']
model_cutoffs['delta_logit_psi_mmsplice_splicemap'] = model_cutoffs['delta_logit_psi_mmsplice']
model_cutoffs['delta_logit_psi_mtsplice_mtsplice'] = model_cutoffs['delta_logit_psi_mmsplice']