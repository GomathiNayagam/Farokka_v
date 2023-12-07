import os
import subprocess
from Bio import SeqIO

def extract_amino_acid_sequences(genbank_file, output_fasta_file):
    """
    Extract amino acid sequences from GenBank file and save in FASTA format.
    :param genbank_file: GenBank file containing query sequences
    :param output_fasta_file: Output FASTA file for amino acid sequences
    :return:
    """
    with open(output_fasta_file, 'w') as fasta_handle, open(genbank_file, 'r') as genbank_handle:
        for record in SeqIO.parse(genbank_handle, 'genbank'):
            # Check if the record has amino acid sequence
            if 'translation' in record.features[0].qualifiers:
                amino_acid_seq = record.features[0].qualifiers['translation'][0]
                seq_record = SeqIO.SeqRecord(Seq(amino_acid_seq), id=record.id, description=record.description)
                SeqIO.write(seq_record, fasta_handle, 'fasta')

def run_diamond_blastp(query_file, database_file, output_file):
    """
    Run diamond blastp search.
    :param query_file: input query file in FASTA format
    :param database_file: Diamond database file
    :param output_file: output file for blastp results
    :return:
    """
    diamond_cmd = f"diamond blastp -d {database_file} -q {query_file} -o {output_file} --evalue 1e-5"

    try:
        subprocess.run(diamond_cmd, shell=True, check=True)
        print("Blastp search completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running Diamond: {e}")

def parse_genbank_and_append_hits(genbank_file, blastp_result_file, output_file):
    """
    Parse query sequences from GenBank file, perform diamond blastp search, and append hit information.
    :param genbank_file: GenBank file containing query sequences
    :param blastp_result_file: blastp result file
    :param output_file: output GenBank file with appended hit information
    :return:
    """
    # Dictionary to store blastp hits
    blastp_hits = {}

    # Read blastp results and store hits in the dictionary
    with open(blastp_result_file, 'r') as blastp_file:
        for line in blastp_file:
            fields = line.strip().split('\t')
            query_id, hit_id = fields[0], fields[1]
            blastp_hits[query_id] = hit_id

    # Read GenBank file, append hit information, and write to output file
    with open(output_file, 'w') as output_handle, open(genbank_file, 'r') as genbank_handle:
        for record in SeqIO.parse(genbank_handle, 'genbank'):
            # Get the GenBank accession or identifier
            genbank_id = record.id

            # Check if the GenBank ID is in the blastp hits dictionary
            if genbank_id in blastp_hits:
                hit_id = blastp_hits[genbank_id]
                record.description += f" (Blastp Hit: {hit_id})"
            SeqIO.write(record, output_handle, 'genbank')

# Example usage:
query_genbank_file_path = "/path/to/your/query_genbank_file.gb"
diamond_db_path = "/path/to/your/diamond_db.dmnd"
blastp_result_path = "/path/to/your/blastp_results.txt"
output_genbank_result_path = "/path/to/your/output_genbank_result.gb"
query_fasta_file_path = "/path/to/your/query.fasta"

# Extract amino acid sequences from GenBank file and save in FASTA format
extract_amino_acid_sequences(query_genbank_file_path, query_fasta_file_path)

# Run diamond blastp search
run_diamond_blastp(query_fasta_file_path, diamond_db_path, blastp_result_path)

# Append hit information to the GenBank file
parse_genbank_and_append_hits(query_genbank_file_path, blastp_result_path, output_genbank_result_path)
