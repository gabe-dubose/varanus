# Function to iteratively annotate variants
def annotate_variants(variants_information, features_database, genome, codon_table):
    import bisect
    import varanus.variant_annotation_utils

    #overall annotations are returned as a dictionary structured as: 
        # annotations = {'Chromosome' : {'Position' : {reference>alternative : [[annotations], [nucleotide_change], [amino_acid_change], [feature_id], [feature_types], [feature_heirarchy], [grand_position]]}}}
    annotations = {}

    for variant in variants_information:

        variant_chromosome = variant[0]
        variant_position = variant[1]
        reference_sequence = variant[2]
        variant_sequence = variant[3]

        #variants attributes are structured as:
        # variant_attributes = [[feature_id(s)], [feature_type(s)], [feature_heirarchy], [annotation_information], [variant_annotation]]
        variant_attributes = [[], [], [], [], []]


        #perform binary search (left) to identify the leftmost position the variant 
        #start can be inserted into the coding sequences locations such that
        #the list maintains the same order. Binary search is performed using 
        #the starting positions of each coding sequence.
        #retreive a list of coding sequences locations
        try:
            coding_features = features_database['features_type_data'][variant_chromosome]['CDS']
            cds_placement = bisect.bisect_left(coding_features, [variant_position, ])
            #check if the stop position of the left adjacent feature is greater than the variant position.
            #If this is the case, the variant is located inside of said feature (in this case, coding sequence).
            if cds_placement > 0:
                feature_stop_position = coding_features[cds_placement-1][1]
                if feature_stop_position > variant_position:
                    feature_id = coding_features[cds_placement-1][2]

                    #add feature id
                    variant_attributes[0].extend([feature_id])
                    #add feature types information
                    variant_attributes[1].extend(['CDS', 'exon', 'mRNA', 'gene'])

                    #get and add feature heirarchy information
                    feature_heirarchy = features_database['features_heirarchy'][variant_chromosome][feature_id]
                    immediate_parent = feature_heirarchy[-1]
                    variant_attributes[2].extend(feature_heirarchy)
                    #in some cases, a CDS can come from multiple different parent. Therefore, jsut grab one parent
                    #and get the rest of the heirarchy (the extended heirarchy will be the same for both parents)
                    if ',' in immediate_parent:
                        first_parent = immediate_parent.split(',')[0]
                        extended_heirarchy = features_database['features_heirarchy'][variant_chromosome][first_parent]
                        variant_attributes[2].extend(extended_heirarchy)

                    #coding sequences that combine to form a protein product should be grouped together
                    #under the same parent, and their order is stored in the "features_order" section
                    #of the features databae. Therefore, the position of the coding sequence can be obtained
                    #by querrying the immediate parent id against the feature order dictionary. 
                    immediate_parent_features = features_database['features_order'][immediate_parent]['CDS']
                    #Find feature ID feature parent (make sure feature ID matches and variant is between start and stop)
                    #This is necessary to annotate start codon variants (i.e., need to location with coding sequence is first)
                    for i in range(len(immediate_parent_features)):
                        if feature_id in immediate_parent_features[i] and variant_position > immediate_parent_features[i][0] and variant_position < immediate_parent_features[i][1]:
                            feature_position_in_parent = [i, len(immediate_parent_features)]
                    #get data needed to annotate the variant
                    feature_start = coding_features[cds_placement-1][0]
                    cds_information = varanus.variant_annotation_utils.get_cds_information(genome, variant_chromosome, variant_position, variant_sequence, reference_sequence, feature_start, codon_table)

                    variant_attributes[3].extend(cds_information)
                    #add cds information

                    #get variant annotaiton
                    cds_annotation = varanus.variant_annotation_utils.get_cds_variant_annotation(variant_attributes, feature_position_in_parent)
                    variant_attributes[4].extend(cds_annotation)

        #skip if sequence doesn't have any CDS
        except:
            pass

        try:
            #if variant is not in a coding sequence, check to see if it is in an exon
            if 'CDS' not in variant_attributes[1]:
                exon_features = features_database['features_type_data'][variant_chromosome]['exon']
                exon_placement = bisect.bisect_left(exon_features, [variant_position, ])
                if exon_placement != 0:
                    feature_stop_position = exon_features[exon_placement-1][1]
                    if feature_stop_position > variant_position:
                        feature_id = exon_features[exon_placement-1][2]

                        #add feature ID, feature types, and feature heirarchy
                        variant_attributes[0].extend([feature_id])
                        variant_attributes[1].extend(['exon', 'mRNA', 'gene'])
                        feature_heirarchy = features_database['features_heirarchy'][variant_chromosome][feature_id]
                        variant_attributes[2].extend(feature_heirarchy)
                    
                        #if variant is not in a coding sequence but is in an exon, variant is likely in the UTR
                        #to confirm, check if variant is in the first exon (5 UTR) or last exon (3 UTR)

                        #get the parent mRNA
                        immediate_parent = feature_heirarchy[-1]
                        #get the order of the exons in the parent mRNA
                        immediate_parent_exons = features_database['features_order'][immediate_parent]['Exons']

                        for exon in immediate_parent_exons:

                            exon_start = exon[0]
                            exon_stop = exon[1]
                            exon_position = immediate_parent_exons.index(exon)
                            
                            #if variant is before first exon
                            if variant_position < exon_start and exon_position == 0:
                                utr_region = '5_prime'
                                variant_attributes[1].extend(['5_UTR'])
                            elif variant_position > exon_stop and exon_position == len(immediate_parent_exons)-1:
                                variant_attributes[1].extend(['3_UTR'])
      
        #skip if exon is not found
        except:
            pass

        try:
            #if variant is not in exon or coding sequence, check to see if it is in mRNA
            if 'exon' not in variant_attributes[1]:
                mRNA_features = features_database['features_type_data'][variant_chromosome]['mRNA']
                mRNA_placement = bisect.bisect_left(mRNA_features, [variant_position, ])
                if mRNA_placement != 0:
                    feature_stop_position = mRNA_features[mRNA_placement-1][1]
                    if feature_stop_position > variant_position:
                        feature_id = mRNA_features[mRNA_placement-1][2]
                        variant_attributes[1].extend(['mRNA', 'gene'])
                        variant_attributes[0].extend([feature_id])
                        feature_heirarchy = features_database['features_heirarchy'][variant_chromosome][feature_id]
                        variant_attributes[2].extend(feature_heirarchy)

                        #if variant is in an mRNA but not an exon, it is likely in an untranslated region.
                        
                        #To identify if variant is in a 5 prime UTR or a 3 prime UTR: see if the variant position is 
                        # before the first exon or after the last exon

                        #query mRNA id against features_order dictionary to get a list of lists in order by feature start
                        mRNA_features_ordered = features_database['features_order'][feature_id]['Exons']
                        #identify if variant position is before the first exon (5_UTR) or after the last exon (3_UTR)
                        for exon_feature in mRNA_features_ordered:
                            exon_start = exon_feature[0]
                            exon_stop = exon_feature[1]
                            exon_position = mRNA_features_ordered.index(exon_feature)
                            #if variant position is before of the first exon
                            if variant_position < exon_start and exon_position == 0:
                                utr_region = 'five_prime'
                                variant_attributes[1].extend(['5_UTR'])
                            #if variant position is after last exon
                            elif variant_position > exon_stop and exon_position == len(mRNA_features_ordered)-1:
                                variant_attributes[1].extend(['3_UTR'])

        except:
            pass

        try:
        #if variant is not in an mRNA, exon, or coding sequence, check to see if it is in a gene
            if 'mRNA' not in variant_attributes[1]:
                gene_features = features_database['features_type_data'][variant_chromosome]['gene']
                gene_placement = bisect.bisect_left(gene_features, [variant_position, ])
                if gene_placement != 0:
                    feature_stop_position = gene_features[gene_placement-1][1]
                    if feature_stop_position > variant_position:
                        variant_attributes[1].extend(['gene'])
                        variant_attributes[0].extend([feature_id])
        except:
            pass

        #if variant is in a CDS, exon, or mRNA, try to find its position in the gene
        #note: Need to adjut to if it is just in a gene at all
        if 'CDS' in variant_attributes[1] or 'exon' in variant_attributes[1] or 'mRNA' in variant_attributes[1]:
            try:
                grand_feature = variant_attributes[2][-1]
                grand_feature_information = features_database['raw_features_data'][variant_chromosome][grand_feature]
                grand_feature_start = grand_feature_information[1]
                grand_feature_stop = grand_feature_information[2]
                variant_grand_position = variant_position - grand_feature_start + 1
            except:
                variant_grand_position = 'NA'
        else:
            variant_grand_position = 'NA'

        try:
        #if variant is not in a gene, variant is intergenic
            if 'gene' not in variant_attributes[1]:
                variant_attributes[1].extend(['intergenic'])
                variant_attributes[0].extend(['No_features'])

                try:
                    gene_features = features_database['features_type_data'][variant_chromosome]['gene']
                    position_in_features = bisect.bisect_left(gene_features, [variant_position, ])

                    #get upstream and downstream genes if feature is intergenic
                    try:
                        upstream_gene = gene_features[position_in_features-1]
                    except:
                        upstream_gene = 'NA'
                    
                    try:
                        downstream_gene = gene_features[position_in_features]
                    except:
                        upstream_gene = ''
                    
                    intergenic_annotation = varanus.variant_annotation_utils.get_intergenic_variant_annotation(upstream_gene, downstream_gene, variant_position)

                    variant_attributes[4].extend([intergenic_annotation])
                except:
                    variant_attributes[4].extend(['intergenic_no_annotation'])
        except:
            pass
            
        #iterate through the remaining feature types to collect more information about 
        #what features the variant is in
        feature_types = list(features_database['features_type_data'][variant_chromosome].keys())
        try:
            feature_types.remove('CDS')
        except:
            pass
        try:
            feature_types.remove('exon')
        except:
            pass
        try:
            feature_types.remove('mRNA')
        except:
            pass
        try:
            feature_types.remove('gene')
        except:
            pass
        try:
            feature_types.remove('region')
        except:
            pass

        for feature_type in feature_types:
            try:
                features = features_database['features_type_data'][variant_chromosome][feature_type]
                feature_placement = bisect.bisect_left(features, [variant_position, ])
                if feature_placement != 0:
                    feature_stop_position = features[feature_placement-1][1]
                    if feature_stop_position > variant_position:
                        feature_id = features[feature_placement-1][2]
                        variant_attributes[1].extend([feature_type])
                        variant_attributes[0].extend([feature_id])
            except:
                pass

        #if variant is not on CDS or Exon, identify if it is in an intron
        if len(variant_attributes[1]) == 2 and 'mRNA' in variant_attributes[1] and 'gene' in variant_attributes[1]:
            variant_attributes[1].extend(['intron'])
            intron_variant_annotation = varanus.variant_annotation_utils.get_intron_variant_annotation()
            variant_attributes[4].extend([intron_variant_annotation])
            

        #set up dictonary
        if variant_chromosome not in annotations:
            annotations[variant_chromosome] = {}
        if str(variant_position) not in annotations[variant_chromosome]:
            annotations[variant_chromosome][str(variant_position)] = {}
        variant_key = f"{reference_sequence}>{variant_sequence}"
        if variant_key not in annotations[variant_chromosome][str(variant_position)]:
            annotations[variant_chromosome][str(variant_position)][variant_key] = []

        #add annotations
        annotations[variant_chromosome][str(variant_position)][variant_key].append(variant_attributes[4])
        #add nucleotide change
        annotations[variant_chromosome][str(variant_position)][variant_key].append([f"{reference_sequence}{variant_position}{variant_sequence}"])
        #add amino acid change if variant is in a coding sequence
        if 'CDS' in variant_attributes[1]:
            #if reference codon and alternative codon are complete codon(s) (sets of 3), add their amino acid changes
            if len(variant_attributes[3][0][0]) % 3 == 0 and len(variant_attributes[3][1][0]) % 3 == 0:
                amino_acid_change = f"{variant_attributes[3][0][1]}{variant_attributes[3][2][3]}{variant_attributes[3][1][1]}"
                annotations[variant_chromosome][str(variant_position)][variant_key].append([amino_acid_change])
        else:
            annotations[variant_chromosome][str(variant_position)][variant_key].append(['NA'])
        
        #add feature id(s)
        annotations[variant_chromosome][str(variant_position)][variant_key].append(variant_attributes[0])
        #add feature type(s)
        annotations[variant_chromosome][str(variant_position)][variant_key].append(variant_attributes[1])
        #add feature heirarchy
        annotations[variant_chromosome][str(variant_position)][variant_key].append(variant_attributes[2])
        #add grand position
        annotations[variant_chromosome][str(variant_position)][variant_key].append(variant_grand_position)


    return annotations