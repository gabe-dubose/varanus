
#function to write database to file
def write_database_file(database, filename):
    import json

    with open(filename, 'w') as outfile:
        json.dump(database, outfile)

#function to write variant annotation to csv file
#Note: annotations = {'Chromosome' : {'Position' : {reference>alternative : [[annotations], [nucleotide_change], [amino_acid_change], [feature_types], [feature_id], [feature_heirarchy]]}}}
#Annotation lines are returned as: chromosme, position, variant, amino acid change, feature id, feature types, feature heirarchy, annotations->
def write_annotations_delimited(variant_annotations, outfile, delimiter):
    with open(outfile, 'a') as outfile:
        for chromosome in variant_annotations:
            for position in variant_annotations[chromosome]:
                for variant in variant_annotations[chromosome][position]:
                    annotation_info = variant_annotations[chromosome][position][variant]
                    
                    #assemble amino acid change field
                    if annotation_info[2][0] == 'NA':
                        amino_acid_change = 'NA:non_coding'
                    else:
                        amino_acid_change = annotation_info[2][0]
                    
                    #assemble feature id field
                    try:
                        if annotation_info[3][0] == 'No_features':
                            try:
                                feature_id = annotation_info[3][1]
                            except:
                                feature_id = 'NA:non_coding'
                        else:
                            feature_id = annotation_info[3][0]
                    except:
                        feature_id = 'NA:non_coding'
                    
                    #assemble feature types field
                    feature_types = ":".join(annotation_info[4])

                    #assemble feature heirarchy field
                    if len(annotation_info[5]) > 0:
                        features_heirarchy = [feature for feature in annotation_info[5] if feature != 'NA']
                        features_heirarchy = ":".join(features_heirarchy)
                    else:
                        features_heirarchy = 'NA:no_heirarchy'

                    #assemble annotation field
                    try:
                        annotation_field = []
                        annotation_s = annotation_info[0][0]
                        if type(annotation_s) == list:
                            for annotation in annotation_s:
                                annotation_to_add = ':'.join(annotation)
                                annotation_field.append(annotation_to_add)
                            annotation_field = ",".join(annotation_field)
                        else:
                            annotation_field = annotation_s
                    except:
                        annotation_field = 'NA:no_additional_annotations'

                    #assmeble and write annotation line
                    annotation_line = f"{chromosome}{delimiter}{position}{delimiter}{variant}{delimiter}{amino_acid_change}{delimiter}{feature_types}{delimiter}{feature_id}{delimiter}{features_heirarchy}{delimiter}{annotation_field}"
                    outfile.write(f"{annotation_line}\n")
