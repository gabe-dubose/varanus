import varanus.database_builder
import varanus.file_loader
import varanus.variant_annotator
import varanus.file_writer
import varanus.defaults

# a single function to perform variant annotation (not usually recommended because this method neglects additionl checks)
def get_variant_annotations(variants_file, features_file, genome_file, codon_table=varanus.defaults.standard_codon_table):
    variants = varanus.file_loader.read_vcf(variants_file)
    database = varanus.database_builder.build_gff3_database(features_file)
    genome = varanus.file_loader.read_genome(genome_file)
    annotations = varanus.variant_annotator.annotate_variants(variants, database, genome, codon_table)
    return annotations

#function to read variants file
def read_variants_file(variants_file):
    variants = varanus.file_loader.read_vcf(variants_file)
    return variants

#function to build features database
def build_database(features_file):
    database = varanus.database_builder.build_gff3_database(features_file) 
    return database

def write_database(database, outfile='varanus_db.json'):
    varanus.file_writer.write_database_file(database, outfile)  

#function to read genome file
def read_genome_file(genome_file):
    genome = varanus.file_loader.read_genome(genome_file)
    return genome

def read_codon_table(codon_table):
    return codon_table

#function to annotate variants
def annotate_variants(variants, database, genome, codon_table = varanus.defaults.standard_codon_table):
    variant_annotations = varanus.variant_annotator.annotate_variants(variants, database, genome, codon_table)
    return variant_annotations
