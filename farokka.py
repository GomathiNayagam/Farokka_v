from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess
import tempfile

def call_orfs(input_sequence, output_fasta):
    # Use PyroDIGAL to predict ORFs
    subprocess.run(["pyrodigal-gv", "-i", input_sequence, "-o", output_fasta])

def translate_orfs(input_fasta, output_translated_fasta):
    # Translate ORFs
    records = list(SeqIO.parse(input_fasta, "fasta"))
    translated_records = [SeqRecord(record.seq.translate(), id=record.id, description="") for record in records]

    # Write translated sequences to a new FASTA file
    SeqIO.write(translated_records, output_translated_fasta, "fasta")

def blastp_search(input_translated_fasta, database_file, output_blast_results):
    # Perform DIAMOND BLASTp search
    subprocess.run(["diamond", "blastp", "-q", input_translated_fasta, "-d", database_file, "--out", output_blast_results])

def create_genbank_file(input_sequence, input_translated_fasta, output_genbank):
    # Load original sequence
    sequence_record = SeqIO.read(input_sequence, "fasta")

    # Load translated sequences
    translated_records = list(SeqIO.parse(input_translated_fasta, "fasta"))

    # Add translated sequences as CDS features in GenBank format
    for i, translated_record in enumerate(translated_records):
        feature_location = f"{i+1}..{i+1+len(translated_record.seq)*3}"  # Assume 1-based coordinates
        feature = f'     CDS             {feature_location}\n                     /translation="{str(translated_record.seq)}"\n'
        sequence_record.seq += feature

    # Write the GenBank file
    SeqIO.write(sequence_record, output_genbank, "genbank")

# Replace 'input_sequence.fasta', 'your_database.fasta', and 'output_prefix' with your actual file names
input_sequence = "input_sequence.fasta"
database_file = "your_database.fasta"
output_prefix = "output_prefix"

# Temporary files
temp_orfs_fasta = tempfile.NamedTemporaryFile(suffix=".fasta").name
temp_translated_fasta = tempfile.NamedTemporaryFile(suffix=".fasta").name
temp_blast_results = tempfile.NamedTemporaryFile(suffix=".tsv").name

# Step 1: Call ORFs
call_orfs(input_sequence, temp_orfs_fasta)

# Step 2: Translate ORFs
translate_orfs(temp_orfs_fasta, temp_translated_fasta)

# Step 3: Perform DIAMOND BLASTp search
blastp_search(temp_translated_fasta, database_file, temp_blast_results)

# Step 4: Create GenBank file
create_genbank_file(input_sequence, temp_translated_fasta, f"{output_prefix}_result.gbk")

# Clean up temporary files
for temp_file in [temp_orfs_fasta, temp_translated_fasta, temp_blast_results]:
    subprocess.run(["rm", temp_file])
