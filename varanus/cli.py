#!/usr/bin/env python3

import argparse
import sys
import datetime
import varanus.api
import varanus.help_messages

parser = argparse.ArgumentParser(add_help=False)
parser._positionals.title = 'varanus commands'
subparser = parser.add_subparsers(dest='command')
parser.add_argument('-h', '--help', required=False, action='store_true')

annotate = subparser.add_parser('annotate-variants', add_help=False)
builddb = subparser.add_parser('build-database')
checkdb = subparser.add_parser('check-database')

annotate._optionals.title = 'arguments'
annotate.add_argument('-h', '--help', required=False, action='store_true')
annotate.add_argument('-v', '--vcf', dest="vcf_file")
annotate.add_argument('-f', '--features', dest="features_file")
annotate.add_argument('-g', '--genome', dest="genome_file")
annotate.add_argument('-c', '--codon_table', dest="codon_table")
annotate.add_argument('-o', '--output', dest="annotations_outfile")
annotate.add_argument('-w', '--writedb', dest="write_database", action="store_true")
annotate.add_argument('-n', '--db_file_name', dest="database_name", default="varanus_db.json")
annotate.add_argument('--verbose', dest="verbose", default="on")

builddb._optionals.title = 'arguments'
builddb.add_argument('-g', '--genome', dest="genome_file")
builddb.add_argument('-f', '--features', dest="features_file")
builddb.add_argument('-o', '--output', dest="database_outfile")

checkdb._optionals.title = 'arguments'
checkdb.add_argument('-d', '--database')
checkdb.add_argument('-o', '--output')

args = parser.parse_args()

#print overall help message
if args.command == None:
    print(varanus.help_messages.main_help_message)

#print help message for annotate-variants command
if args.command == 'annotate-variants' and len(sys.argv) == 2:
    print(varanus.help_messages.annotate_variants_help_message)

if args.command == 'annotate-variants' and args.help == True:
    print(varanus.help_messages.annotate_variants_help_message)

#Command line interface for annotating variants

#annotate-variants function
def annotate_variants():

    #read variants file
    #verbose options
    if args.verbose == 'on':
        print("Reading variants file...")
        variant_reading_start = datetime.datetime.now()
    #read vcf file
    variants = varanus.api.read_variants_file(args.vcf_file)
    if args.verbose == 'on':
        print(f"Variants read in {datetime.datetime.now() - variant_reading_start}")

            
    #building features database
    #verbose options
    if args.verbose == 'on':
        print("Building database...")
        db_build_start = datetime.datetime.now()
    #build database
    database = varanus.api.build_database(features_file=args.features_file)
    #verbose options
    if args.verbose == 'on':
        print(f"Database built in {datetime.datetime.now() - db_build_start}")
    
    # (optional) writing features database
    if args.write_database == True or args.database_name != 'varanus_db.json':
        #verbose options
        if args.verbose == 'on':
            print("Writing database...")
            writing_database_start = datetime.datetime.now()
        #write database to file
        varanus.api.write_database(database, args.database_name)
        #verbose options
        if args.verbose == 'on':
            print(f"Database written in {datetime.datetime.now() - writing_database_start}")
    
    #read genome file
    #verbose option
    if args.verbose == 'on':
        print("Reading genome...")
        genome_reading_start = datetime.datetime.now()
    #read genome
    genome = varanus.api.read_genome_file(args.genome_file)
    #verbose option
    if args.verbose == 'on':
        print(f"Genome read in {datetime.datetime.now() - genome_reading_start}")

    #change codon table if specificed
    if bool(args.codon_table) == True:
        codon_table = varanus.api.read_codon_table(args.codon_table)
    #verbose option
        if args.verbose == True:
            print("Codon table loaded.")

    #annotate variants
    #verbose option
    if args.verbose == 'on':
        print("Annotating variants...")
        annotation_start = datetime.datetime.now()
    if bool(args.codon_table) == True:
        variant_annotations = varanus.api.annotate_variants(variants, database, genome, codon_table=codon_table)
    else:
        variant_annotations = varanus.api.annotate_variants(variants, database, genome)
    #verbose option
    if args.verbose == 'on':
        print(f"Variants annotated in {datetime.datetime.now() - annotation_start}")

#annotate-variants
if args.command == 'annotate-variants' and len(sys.argv) > 2:
    annotate_variants()