# This module contains default parameters used throughout varanus

#default feature_id regex terms
feature_id_regex_terms = ['ID=.*?;', 'ID=', ';']

#default parent_id regex terms
parent_id_regex_terms = ['Parent=.*?;|Parent=.*?\n', 'Parent=', ';', '\n']

#default gene name regex terms
gene_name_regex_terms = ['gene=.*?;', 'gene=', ';']

#default protein product regex terms
protein_product_regex_terms = ['product=.*?;', 'product=', ';']

#default locus tag regex terms
locus_tag_regex_terms = ['locus_tag=.*?;', 'locus_tag=', ';']

#Total default regex terms
default_info_search_terms = {'feature_id_terms' : feature_id_regex_terms, 
    'parent_id_terms' : parent_id_regex_terms, 'gene_name_terms' : gene_name_regex_terms,
    'protein_product_terms' : protein_product_regex_terms, 'locus_tag_terms' : locus_tag_regex_terms}

#standard codon table structured as:
#    standard_codon_table = {'CODON' : ['AMINO_ACID', 'BIOCHEMISTRY', 'INFO']...}
standard_codon_table = {'TTT': ['PHE', 'AROMATIC_R_GROUP', 'NA'], 'TTC': ['PHE', 'AROMATIC_R_GROUP', 'NA'], 
    'TTA': ['LEU', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 'TTG': ['LEU', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 
    'TCT': ['SER', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 'TCC': ['SER', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 
    'TCA': ['SER', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 'TCG': ['SER', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 
    'TAT': ['TYR', 'AROMATIC_R_GROUP', 'NA'], 'TAC': ['TYR', 'AROMATIC_R_GROUP', 'NA'], 
    'TAA': ['NA', 'NA', 'STOP'], 'TAG': ['NA', 'NA', 'STOP'], 'TGT': ['CYS', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 
    'TGC': ['CYS', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 'TGA': ['NA', 'NA', 'STOP'], 'TGG': ['TRP', 'AROMATIC_R_GROUP', 'NA'], 
    'CTT': ['LEU', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 'CTC': ['LEU', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 
    'CTA': ['LEU', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 'CTG': ['LEU', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 
    'CCT': ['PRO', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 'CCC': ['PRO', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 
    'CCA': ['PRO', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 'CCG': ['PRO', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 
    'CAT': ['HIS', 'POSITIVE_CHARGED_R_GROUP', 'NA'], 'CAC': ['HIS', 'POSITIVE_CHARGED_R_GROUP', 'NA'], 
    'CAA': ['GLN', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 'CAG': ['GLN', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 
    'CGT': ['ARG', 'POSITIVE_CHARGED_R_GROUP', 'NA'], 'CGC': ['ARG', 'POSITIVE_CHARGED_R_GROUP', 'NA'], 
    'CGA': ['ARG', 'POSITIVE_CHARGED_R_GROUP', 'NA'], 'CGG': ['ARG', 'POSITIVE_CHARGED_R_GROUP', 'NA'], 
    'ATT': ['ILE', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 'ATC': ['ILE', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 
    'ATA': ['ILE', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 'ATG': ['MET', 'NONPOLAR_ALIPHATIC_R_GROUP', 'START'], 
    'ACT': ['THR', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 'ACC': ['THR', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 'ACA': ['THR', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 
    'ACG': ['THR', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 'AAT': ['ASN', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 'AAC': ['ASN', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 
    'AAA': ['LYS', 'POSITIVE_CHARGED_R_GROUP', 'NA'], 'AAG': ['LYS', 'POSITIVE_CHARGED_R_GROUP', 'NA'], 'AGT': ['SER', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 
    'AGC': ['SER', 'POLAR_UNCHARGED_R_GROUP', 'NA'], 'AGA': ['ARG', 'POSITIVE_CHARGED_R_GROUP', 'NA'], 'AGG': ['ARG', 'POSITIVE_CHARGED_R_GROUP', 'NA'], 
    'GTT': ['VAL', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 'GTC': ['VAL', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 'GTA': ['VAL', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 
    'GTG': ['VAL', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 'GCT': ['ALA', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 'GCC': ['ALA', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 
    'GCA': ['ALA', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 'GCG': ['ALA', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 'GAT': ['ASP', 'NEGATIVE_CHARGED_R_GROUP', 'NA'], 
    'GAC': ['ASP', 'NEGATIVE_CHARGED_R_GROUP', 'NA'], 'GAA': ['GLU', 'NEGATIVE_CHARGED_R_GROUP', 'NA'], 'GAG': ['GLU', 'NEGATIVE_CHARGED_R_GROUP', 'NA'], 
    'GGT': ['GLY', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 'GGC': ['GLY', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 'GGA': ['GLY', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA'], 
    'GGG': ['GLY', 'NONPOLAR_ALIPHATIC_R_GROUP', 'NA']}