# Function to get the information necessary to annotate coding sequence variant
# Information needed to annotate variants includes:
# 1) Reference information: Reference codon, reference amino acid, reference biochemistry, 
#   if reference is start or stop codon (True or False), and position in codon (to assemble alternative codon)
# 2) Alternative information: Alternative codon, alternative amino acid, alternative biochemistry,
#   if alternative is start or stop codon (True or False), 
# 3) Positional information: Previous codon start or stop, next codon start or stop, position in gene
def get_cds_information(genome, chromosome, variant_position, variant_sequence, reference_sequence, feature_start, codon_table):
    import math

    warning_messages = []

    position_in_feature = variant_position - feature_start + 1
    codon_position = math.ceil(position_in_feature/3)
    position_in_codon = position_in_feature%3 + 1

    #get reference codon based on the position of the variant in its chromosome
    if position_in_codon == 1:
        reference_codon = str(genome[chromosome][variant_position-1:variant_position+2]).upper()
        #check to make sure reference from VCF matches reference found in genome file
        if reference_codon[0] != reference_sequence:
            warning_messages.append(f"Warning: Reference sequence from VCF '{reference_sequence}' does not match reference sequence found in genome '{reference_codon[0]}'")
    elif position_in_codon == 2:
        reference_codon = str(genome[chromosome][variant_position-2:variant_position+1]).upper()
        #check to make sure reference from VCF matches reference found in genome file
        if reference_codon[1] != reference_sequence:
            warning_messages.append(f"Warning: Reference sequence from VCF '{reference_sequence}' does not match reference sequence found in genome '{reference_codon[1]}'")
    elif position_in_codon == 3:
        reference_codon = str(genome[chromosome][variant_position-3:variant_position]).upper()
        #check to make sure reference from VCF matches reference found in genome file
        if reference_codon[2] != reference_sequence:
            warning_messages.append(f"Warning: Reference sequence from VCF '{reference_sequence}' does not match reference sequence found in genome '{reference_codon[2]}'")

    #get reference amino acid, biochemistry, and start/stop informaiton
    if reference_codon in codon_table:
        reference_amino_acid = codon_table[reference_codon][0]
        reference_biochemistry = codon_table[reference_codon][1]
        reference_start_stop_bool = codon_table[reference_codon][2]
    else:
        warning_messages.append(f"Error: Reference codon '{reference_codon}' is not a recognized codon.")
        reference_codon = 'Fail'
        reference_amino_acid = 'Fail'
        reference_biochemistry = 'Fail'
        reference_start_stop_bool = 'Fail'

    #assemble reference information
    reference_information = [reference_codon, reference_amino_acid, reference_biochemistry, reference_start_stop_bool]

    #get alternative codon, amino acid, biochemistry, and start/stop information
    replacement_position = position_in_codon - 1
    alternative_codon = reference_codon[:replacement_position] + variant_sequence +  reference_codon[replacement_position+1:]

    #if variant is not insertion or deletion, get alternative information
    if len(alternative_codon) == len(reference_codon):
        if alternative_codon in codon_table:
            alternative_amino_acid = codon_table[alternative_codon][0]
            alternative_biochemistry = codon_table[alternative_codon][1]
            alternative_start_stop_bool = codon_table[alternative_codon][2]
        else:
            warning_messages.append(f"Error: Alternative codon '{alternative_codon}' is not a recognized codon.")
            alternative_codon = 'Fail'
            alternative_amino_acid = 'Fail'
            alternative_biochemistry = 'Fail'
            alternative_start_stop_bool = 'Fail'
    else:
        alternative_amino_acid = 'NA'
        alternative_biochemistry = 'NA'
        alternative_start_stop_bool = 'NA'

    #assemble alternative information
    alternative_information = [alternative_codon, alternative_amino_acid, alternative_biochemistry, alternative_start_stop_bool]

    #try to get positional information
    #get previous codon information
    try:
        if codon_position > 1:
            if position_in_codon == 1:
                previous_codon = str(genome[chromosome][variant_position-4:variant_position-1]).upper()
            elif position_in_codon == 2:
                previous_codon = str(genome[chromosome][variant_position-5:variant_position-2]).upper()
            elif position_in_codon == 3:
                previous_codon = str(genome[chromosome][variant_position-6:variant_position-3]).upper()
            #make sure codon in is a valid codon
            if previous_codon in codon_table and len(previous_codon) == 3:
                previous_start_stop_bool = codon_table[previous_codon][2]
            else:
                warning_messages.append(f"Warning: Previous codon '{previous_codon}' is not a recognized codon.")
                previous_start_stop_bool = 'Fail'
        else:
            previous_start_stop_bool = 'NA'
    except:
        previous_start_stop_bool = 'NA'
    
    #get next codon information
    try:
        if position_in_codon == 1:
            next_codon = str(genome[chromosome][variant_position+2:variant_position+5]).upper()
        elif position_in_codon == 2:
            next_codon = str(genome[chromosome][variant_position+3:variant_position+6]).upper()
        elif position_in_codon == 3:
            next_codon = str(genome[chromosome][variant_position+4:variant_position+7]).upper()
        #make sure codon in is a valid codon
        if next_codon in codon_table and len(next_codon) == 3:
            next_start_stop_bool = codon_table[next_codon][2]
        else:
            next_start_stop_bool = 'Fail'
    except:
        next_start_stop_bool = 'NA'
    
    #assemble positional information
    positional_information = [previous_start_stop_bool, next_start_stop_bool, position_in_feature, codon_position, position_in_codon]

    #assemble return information
    cds_information = [reference_information, alternative_information, positional_information, warning_messages]

    return cds_information

#Function to get coding sequence variant annotation
def get_cds_variant_annotation(variant_attributes, feature_position_in_parent):

    variant_annotation = []

    #unpack information for readability
    cds_information = variant_attributes[3]
    reference_codon = cds_information[0][0]
    reference_amino_acid = cds_information[0][1]
    reference_biochemistry = cds_information[0][2]
    reference_start_stop_bool = cds_information[0][3]

    alternative_codon = cds_information[1][0]
    alternative_amino_acid = cds_information[1][1]
    alternative_biochemistry = cds_information[1][2]
    alternative_start_stop_bool = cds_information[1][3]

    previous_start_stop_bool = cds_information[2][0]
    next_start_stop_bool = cds_information[2][1]
    position_in_feature = cds_information[2][2]
    codon_position = cds_information[2][3]
    position_in_codon = cds_information[2][4]

    warning_messages = cds_information[3]


    #synonymous annotation parent
    if reference_amino_acid == alternative_amino_acid:
        #stop_retained annotation
        if reference_start_stop_bool == 'STOP' and alternative_start_stop_bool == 'STOP':
            variant_annotation.append(['synonymous', 'stop_retained'])
        #synonymous annotation
        else:
            variant_annotation.append(['synonymous'])
            
    #nonsynonymous annotation aprent
    elif reference_amino_acid != alternative_amino_acid:
        #missense annotation parent
        if len(alternative_codon) == len(reference_codon):
            #conservative_missense annotation
            if reference_biochemistry == alternative_biochemistry:
                variant_annotation.append(['nonsynonymous','missense','conservative_missense_variant'])
            #nonconservative_missense_variant annotation
            elif reference_biochemistry != alternative_biochemistry:
                variant_annotation.append(['nonsynonymous','missense','nonconservative_missense_variant'])
            #redundant_inserted_stop_gained annotation
            elif alternative_start_stop_bool == 'STOP' and previous_start_stop_bool == 'STOP' or alternative_start_stop_bool == 'STOP' and next_start_stop_bool == 'STOP':
                variant_annotation.append(['nonsynonymous','missense','redundant_inserted_stop_gained'])
                #start_lost annotation
            elif reference_start_stop_bool == 'START' and alternative_start_stop_bool != 'START' and position_in_feature == 1 and feature_position_in_parent[0] == 0:
                variant_annotation.append(['nonsynonymous','missense','start_lost'])
            else:
                #missense annotation
                variant_annotation.append(['nonsynonymous','missense'])
        #inframe_deletion annotation parent
        if len(alternative_codon) < len(reference_codon):
            #conservative_inframe_deletion annotation
            if len(reference_codon) - len(alternative_codon) % 3 == 0 and position_in_codon == 1:
                variant_annotation.append(['nonsynonymous', 'inframe_deletion', 'conservative_inframe_deletion'])
            #disruptive_inframe_deletion annotation
            if len(reference_codon) - len(alternative_codon) % 3 != 0 and position_in_codon != 1:
                variant_annotation.append(['nonsynonymous', 'inframe_deletion', 'disruptive_inframe_deletion'])
            #inframe_deletion annotation
            else:
                variant_annotation.append(['nonsynonymous', 'inframe_deletion'])

        #inframe_insertion annotation parent
        if len(alternative_codon) > len(reference_codon):
            #conservative_inframe_insertion annotation
            if len(alternative_codon) - len(reference_codon) % 3 == 0 and position_in_codon == 1 or position_in_codon == 3:
                variant_annotation.append(['nonsynonymous', 'inframe_insertion', 'conservative_inframe_insertion'])
            #disruptive_inframe_insertion
            elif len(alternative_codon) - len(reference_codon) % 3 == 0 and position_in_codon == 2:
                variant_annotation.append(['nonsynonymous', 'inframe_insertion', 'disruptive_inframe_insertion'])
            #inframe_insertion annotation
            else:
                variant_annotation.append(['nonsynonymous', 'inframe_insertion'])

        #nonsynonymous annotation
        if len(variant_annotation) == 0:
            variant_annotation.append(['nonsynonymous'])

    annotation = [variant_annotation, warning_messages]
    return annotation


#Function to get UTR variant annotation
def get_utr_variant_annotation():
    pass

#Function to get intron variant annotation
def get_intron_variant_annotation():
    pass

#Function to get intergenic variant annotation
def get_intergenic_variant_annotation():
    pass