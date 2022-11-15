main_help_message = '''usage: varanus [command] [options]

commands:
annotate-variants   
package-database
check-database

optional arguments: 
-h, --help  Show help message and exit
'''

annotate_variants_help_message = '''usage: varanus annotate-variants [options]
minimum requirements:
option 1:
-v, --vcf\t\t\tvariant call file in VCF format
-f, --features\t\t\tfeatures file in GFF or GTF format
-g, --genome\t\t\tgenome file in FASTA format

option 2:
-v, --vcf\t\t\tvariant call file in VCF format
--varanus-database-artifact\tvaranus database generated from the 'package-database' command
'''