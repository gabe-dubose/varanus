#!/usr/bin/env python3

import argparse
import sys
import varanus.help_messages

parser = argparse.ArgumentParser(add_help=False)
parser._positionals.title = 'varanus commands'
subparser = parser.add_subparsers(dest='command')
parser.add_argument('-h', '--help', required=False, action='store_true')

annotate = subparser.add_parser('annotate-variants', add_help=False)
packagedb = subparser.add_parser('package-database')
checkdb = subparser.add_parser('check-database')

annotate._optionals.title = 'arguments'
annotate.add_argument('-h', '--help', required=False, action='store_true')
annotate.add_argument('-v', '--vcf')
annotate.add_argument('-f', '--features')
annotate.add_argument('-g', '--genome')

packagedb._optionals.title = 'arguments'
packagedb.add_argument('-g', '--genome')
packagedb.add_argument('-f', '--features')
packagedb.add_argument('-o', '--output')

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