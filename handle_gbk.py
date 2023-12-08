from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess

def perform_diamond_search(query_sequence, database_file):
    # Write the query sequence to a temporary FASTA file
    temp_fasta_file = "temp_query.fasta"
    SeqIO.write(SeqRecord(query_sequence), temp_fasta_file, "fasta")

    # Perform DIAMOND search
    diamond_result = subprocess.run(
        ["diamond", "blastp", "-q", temp_fasta_file, "-d", database_file, "--outfmt", "6 sseqid"],
        capture_output=True, text=True
    )

    # Remove temporary FASTA file
    subprocess.run(["rm", temp_fasta_file])

    # Parse DIAMOND results
    blast_hits = [line.strip() for line in diamond_result.stdout.split("\n") if line.strip()]

    return blast_hits

def update_genbank_with_blast_hits(genbank_file, blast_hits):
    # Load GenBank file
    records = list(SeqIO.parse(genbank_file, "genbank"))

    # Update the product record in each SeqFeature with the corresponding DIAMOND hit
    for i, record in enumerate(records):
        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                translation = feature.qualifiers["translation"][0]
                if translation in blast_hits:
                    feature.qualifiers["product"] = [blast_hits[translation]]

    return records

# Replace 'input.gbk', 'your_database.fasta', and 'output.gbk' with your actual file names
input_genbank_file = "input.gbk"
database_file = "your_database.fasta"
output_genbank_file = "output.gbk"

# Extract translated sequences from GenBank file
translated_sequences = [feature.qualifiers["translation"][0] for record in SeqIO.parse(input_genbank_file, "genbank") for feature in record.features if feature.type == "CDS" and "translation" in feature.qualifiers]

# Perform DIAMOND search for each translated sequence
blast_hits_dict = {}
for sequence in translated_sequences:
    hits = perform_diamond_search(sequence, database_file)
    if hits:
        blast_hits_dict[sequence] = hits[0]  # Assume the first hit is the most relevant

# Update GenBank file with DIAMOND hits
updated_records = update_genbank_with_blast_hits(input_genbank_file, blast_hits_dict)

# Write the updated GenBank file
SeqIO.write(updated_records, output_genbank_file, "genbank")

print(f"Updated GenBank file '{output_genbank_file}' created successfully.")
