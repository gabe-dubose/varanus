import varanus.defaults
import varanus.database_build_utils
import bisect
#A function to build a features database from a GFF3 file to be used for variant annotation 
# This function is specifically for gff3 files because they include parent-child relationships
# between features, which is necessary to build the features heirarchy. This heirarchy contains 
# parent-child relations between the features in a gff3 file. For this function to
# work correctly, it must be able to locate the feature ID and parent ID for each feature.
# Default regular expression terms are recorded in database_building_utils.
# Each unique feature is recorded as a key in the heirarchy dictionary with a list value 
# of its total lineage (lineage being parents, grandparents, etc.):
#
# This database consists of four primary components:
#
# 1) A dictionary object to hold raw features data structured as:
#   raw_features_data = {'Chromosome' : {'feature_id' : [feature_type, feature_start, feature_stop, feature_strand, feature_phase],...},...}
#   *This allows for information about each feature to be rereived based on the querry of feature identifiers
#
# 2) A dictionary object to hold information about features heirarchy (parent-child relationships), structured as:
#   features_heirarchy = {'Chromosome' : {feature1 : [], feature2 : [feature1], feature3 = [feature1, feature2],...},...}
#   *This allows for each features lineage (parents of parents, etc.) to be querried
#
# 3) A dictionary object to hold ordered information of each feature type, structured as:
#   features_type_data = {'Chromosome' : {'feature_type' : [feature_start, feature_stop, feature_id], [feature_start, feature_stop, feature_id],...,}} 
#   *This is the primary dictionary used to find which features a variant is located in
#
# 4) A dictionary object to hold the order of each coding sequence and exon in each mRNA/gene parent, structured as:
#   features_order = {'parent_id' : {'CDS' : [start, stop, feature_id],...}, 'Exons' : [start, stop, feature_id],...}
#   *This allows for sequential assemblies of each coding sequence or exon to be retreived based on parent_id querries
#
#Note:
# The features heirarchy and features data structures are kept separate because the primary 
# search algoritm implemented to locate what features variants are present in is a binary search;
# therefore, the features locations must be sorted. To my limited knowledge, if feature location 
# information is included in the heirarchy, a tree traversal algorithm would be required, which would 
# negatively impact performance. Also, some of these data structures can most certainly be combined, 
# which I might get back to later. 

def build_gff3_database(features_file):

    raw_features_data = {}
    features_heirarchy = {}
    features_type_data = {}
    features_order = {}

    with open(features_file, 'r') as infile:
        lines = infile.readlines()
    
    for line in lines:
        if line[0] != '#':

            chromosome = line.split('\t')[0]
            feature_type = line.split('\t')[2]
            feature_start = int(line.split('\t')[3])
            feature_stop = int(line.split('\t')[4])
            feature_strand = line.split('\t')[6]
            feature_phase = line.split('\t')[7]
            info_field = line.split('\t')[8]

            feature_id = varanus.database_build_utils.get_feature_id(info_field, feature_type, feature_start, feature_stop)

            #add to raw data
            varanus.database_build_utils.add_raw_data(raw_features_data, chromosome, feature_id, feature_type, feature_start, feature_stop, feature_strand, feature_phase, info_field)

            #add to features heirarchy
            parent_id = varanus.database_build_utils.get_parent_id(info_field, varanus.defaults.parent_id_regex_terms)
            
            #add information to feature heirarchy
            varanus.database_build_utils.add_feature_heirarchy(features_heirarchy, chromosome, feature_id, parent_id)

            #add information to features data dictionary
            varanus.database_build_utils.add_feature_type_data(features_type_data, chromosome, feature_type, feature_start, feature_stop, feature_id)
            #add information to features order
            if feature_type == 'CDS' or feature_type == 'exon':
                if parent_id not in features_order:
                    features_order[parent_id] = {'CDS' : [], 'Exons' : []}
                if feature_type == 'CDS':
                    bisect.insort_left(features_order[parent_id]['CDS'], [feature_start, feature_stop, feature_id])
                elif feature_type == 'exon':
                    bisect.insort_left(features_order[parent_id]['Exons'], [feature_start, feature_stop, feature_id])

    #sort features by start location
    for chromosome in features_type_data:
        for feature_type in features_type_data[chromosome]:
            try:
                features_type_data[chromosome][feature_type].sort()
            except:
                pass

    #return final database
    database = {'raw_features_data' : raw_features_data, 'features_heirarchy' : features_heirarchy, 'features_type_data' : features_type_data, 'features_order' : features_order}
    return database