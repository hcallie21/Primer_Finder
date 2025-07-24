from filecmp import cmp
import gffutils
#sequence I/O
from Bio import SeqIO
from Bio.Seq import Seq
import primer_eval

#PURPOSE: given a gff3 & fasta file, find a suitable primer for a gene of interest
    # Not all fasta files are annotated, should be able to accomodate
gff3 = "./sequence.gff3"
db_file = "./sequence.db"
fastafile = "./sequence_chrI.fasta"

#only for searching for gene name! change to other stuff later?
name_input = input("Enter gene name to search: ").strip()

#create a  SQL database from the GFF3 file
gffutils.create_db(
    gff3,
    dbfn=db_file,
    force=True,
    keep_order=True,
    merge_strategy="merge",
    sort_attribute_values=True,
)

db = gffutils.FeatureDB(db_file, keep_order=True)

output = None
#finds coding sequence index in gff3 file for gene name
for gene in db.features_of_type("CDS"):
    if name_input.lower() in gene.attributes.get("gene", [""])[0].lower():
        gene_name = gene.attributes.get("gene", ["N/A"])[0]
        gene_id = gene.attributes.get("Dbxref", ["N/A"])[0]
        gene_id = gene_id.split(":")[-1]  # Extract the actual ID if it's in a format like "gene:ID"
        output = [gene_id, gene.start, gene.end, gene.strand, gene.attributes]
        break

#return the gene name information   
if output is not None: 
    print(f"Found Gene Name: {name_input}") 
    print(f"Gene ID: {output[0]}, Start: {output[1]}, End: {output[2]}, Strand: {output[3]}")    
else:
    print(f"No gene found with the name: {name_input}")
    exit(1)

#read the fasta file (whole chr) and get the length of the genome
record = SeqIO.read(fastafile, "fasta")
genome_length = len(record)

#if gene is on the minus strand, convert coordinates to plus strand
#CDS_sequence is the sequence of the gene trying to find primers for
if output[3] == "-":
    cds_sequence = record[output[1] - 1 : output[2]].reverse_complement()
else: 
    cds_sequence = record[output[1] - 1 : output[2]]

#snag sequences upstream and downstream of the CDS to use for primer design
region_length = 300  # Length of the region to consider for primer design
upstream_candidate_region = record[output[1] - 30 - region_length: output[1] - 30] # Upstream region
downstream_candidate_region = record[output[2] + 30 : output[2] + 30 + region_length].reverse_complement()

#find suitable primers for the upstream and downstream candidate regions
#min_length and max_length are the minimum and maximum lengths of the primer
up_primer_seq, down_primer_seq = None, None

up_result = primer_eval.find_suitable_primer(upstream_candidate_region, min_length=18, max_length=35)
if up_result is not None:
    up_primer_seq, x = up_result
    # proceed with using the primer
else:
    print("No suitable primer found.")

dwn_result = primer_eval.find_suitable_primer(downstream_candidate_region, min_length=18, max_length=25)
if dwn_result is not None:
    down_primer_seq, y = dwn_result
    # proceed with using the primer
else:
    print("No suitable primer found.")

#make found some primers
if up_primer_seq is not None:
    print(f"Upstream Primer Sequence: {up_primer_seq} at starting index {output[1] -30-region_length + x} on the plus strand of length {len(up_primer_seq)}")
else:
    print("No suitable upstream primer found.")
if down_primer_seq is not None:
    print(f"Downstream Primer Sequence: {down_primer_seq} at starting index {output[2] + region_length + 30 - y} on the plus strand of length {len(down_primer_seq)}")
else: 
    print("No suitable downstream primer found.")

exit(1)
