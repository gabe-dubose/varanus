
# Function to read VCF file and return variants structured as:
#     variants = [[Chromosome, Position, Reference, Alternative],...]
def read_vcf(vcf_file):
    import re

    variants = []
    #read vcf file and extract necessary fields
    with open(vcf_file, 'r') as infile:
        lines = infile.readlines()

    for line in lines:
        comment = bool(re.match('^#', line))
        if comment == False:
            chromosome = line.split('\t')[0]
            position = int(line.split('\t')[1])
            reference = line.split('\t')[3]
            alternative = line.split('\t')[4]

            variants.append([chromosome, position, reference, alternative])
    return variants

# Function to read genome and return it as a pyfaidx fasta object
def read_genome(genome_file):
    import pyfaidx
    genome = pyfaidx.Fasta(genome_file)
    return genome

# Function to read database file in json format
def read_json_database(database_file):
    import json
    with open(database_file, 'r') as infile:
        database = json.load(infile)
    return database

# Function to read custom codon table
def read_codon_table(codon_table):
    pass