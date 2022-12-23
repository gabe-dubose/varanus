
import re
import varanus.defaults
#Function to find the feature ID from a given info field
# Function takes an info field (string) and a list of regular expressions formatted as:
#   regular_expressions = [primary_search_term, trim_term1, trim_term2,...]
# The first regular expression in the list is used as the primary search term 
# to reliably extract the feature ID from the info field. In some instances, additions
# characters will need to be provided to provide a unique enough pattern to reliably perform
# this pattern mathcing. Therefore, the remaining expressions in the list are used to trim excess
# characters that were used for pattern mathching, but not actually part of the feature ID. 
def get_feature_id(info_field, feature_type, feature_start, feature_stop, regular_expressions = varanus.defaults.default_info_search_terms['feature_id_terms']):

    feature_id = re.search(regular_expressions[0], info_field)
    if bool(feature_id) == True:
        feature_id = feature_id.group()
        for trim_term in regular_expressions[1:]:
            feature_id = re.sub(trim_term, '', feature_id)
    #if feature ID is not found, assemble feature ID from feature type, start, and stop
    else:
        feature_id = f"{feature_type}_{feature_start}_{feature_stop}"
    return feature_id

#Function to find the parent ID from a given info field
# This function works the same way as the 'get_feature_id' function
def get_parent_id(info_field, regular_expressions = varanus.defaults.default_info_search_terms['parent_id_terms']):
    
    parent_id = re.search(regular_expressions[0], info_field)
    if bool(parent_id) == True:
        parent_id = parent_id.group()
        for trim_term in regular_expressions[1:]:
            if trim_term in parent_id:
                parent_id = re.sub(trim_term, '', parent_id)
    else:
        parent_id = 'NA'
    return parent_id

#function to get gene name
def get_gene_name(info_field, regular_expressions = varanus.defaults.default_info_search_terms['gene_name_terms']):
    gene_name = re.search(regular_expressions[0], info_field)
    if bool(gene_name) == True:
        gene_name = gene_name.group()
        for trim_term in regular_expressions[1:]:
            gene_name = re.sub(trim_term, '', gene_name)
    else:
        gene_name = 'NA'
    return gene_name

#function to get protein product
def get_protein_product(info_field, regular_expressions = varanus.defaults.default_info_search_terms['protein_product_terms']):
    protein_product = re.search(regular_expressions[0], info_field)
    if bool(protein_product) == True:
        protein_product = protein_product.group()
        for trim_term in regular_expressions[1:]:
            protein_product = re.sub(trim_term, '', protein_product)
    else:
        protein_product = 'NA'
    return protein_product

#function to get locus tag
def get_locus_tag(info_field, regular_expressions = varanus.defaults.default_info_search_terms['locus_tag_terms']):
    locus_tag = re.search(regular_expressions[0], info_field)
    if bool(locus_tag) == True:
        locus_tag = locus_tag.group()
        for trim_term in regular_expressions[1:]:
            locus_tag = re.sub(trim_term, '', locus_tag)
    else:
        locus_tag = 'NA'
    return locus_tag  

#function to add raw feature data dictionary
def add_raw_data(raw_features_data_dict, chromosome, feature_id, feature_type, feature_start, feature_stop, feature_strand, feature_phase, info_field):

    #attempt to get additional info fields
    gene_name = get_gene_name(info_field)
    protein_product = get_protein_product(info_field)
    locus_tag = get_locus_tag(info_field)

    if chromosome not in raw_features_data_dict:
        raw_features_data_dict[chromosome] = {}
    if feature_id not in raw_features_data_dict[chromosome]:
        raw_features_data_dict[chromosome][feature_id] = [feature_type, feature_start, feature_stop, feature_strand, feature_phase, [gene_name, protein_product, locus_tag]]

#Function to add feature information to feature type data dictionary.
# Dictionary is structured as:
# features_type_data = {'Chromosome' : {'feature_type' : [feature_start, feature_stop, feature_id, feature_strand, feature_phase]
#                                                     [feature_start, feature_stop, feature_id, feature_strand, feature_phase],...,}}
def add_feature_type_data(features_type_data, chromosome, feature_type, feature_start, feature_stop, feature_id):
    #add new chromosome
    if chromosome not in features_type_data:
        features_type_data[chromosome] = {}
    #add new feature
    if feature_type not in features_type_data[chromosome]:
        features_type_data[chromosome][feature_type] = []
    #add information to associated features
    features_type_data[chromosome][feature_type].append([feature_start, feature_stop, feature_id])

#Function to add feature to feature heirarchy
def add_feature_heirarchy(features_heirarchy_dict, chromosome, feature_id, parent_id):

    #add chromosome to heirarchy
    if chromosome not in features_heirarchy_dict:
        features_heirarchy_dict[chromosome] = {}

    #if feature is not in dictionary and it doesn't have a parent, add with empty list
    if feature_id not in features_heirarchy_dict[chromosome] and parent_id == None:
        features_heirarchy_dict[chromosome][feature_id] = []
    #if parent ID is in dictionary, add feature with parents lineage
    elif parent_id in features_heirarchy_dict[chromosome]:
        features_heirarchy_dict[chromosome][feature_id] = features_heirarchy_dict[chromosome][parent_id] + [parent_id]
    #if parent ID is not none and is not in dictionary, add it. Then add feature with parents lineage as value
    elif feature_id not in features_heirarchy_dict[chromosome] and parent_id != None:
        features_heirarchy_dict[chromosome][parent_id] = []
        features_heirarchy_dict[chromosome][feature_id] = features_heirarchy_dict[chromosome][parent_id] + [parent_id]
    #note: features with duplicate feature_ids are only added once: