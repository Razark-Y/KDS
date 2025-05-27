import os
import tempfile
from supabase import create_client, Client
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

# ---- CONFIG ----
SUPABASE_URL = "https://unwqzwvqjdcwmlwafxoo.supabase.co"
SUPABASE_KEY = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InVud3F6d3ZxamRjd21sd2FmeG9vIiwicm9sZSI6ImFub24iLCJpYXQiOjE3NDY4NjU0NzUsImV4cCI6MjA2MjQ0MTQ3NX0.7P3a3cL3JPUidguddJo7HFgy25Ft60UKPpfo4aznqlo"
TABLE_NAME = "Nucleotide"
FASTA_COLUMN = "FASTA"
MAFFT_COLUMN = "MAFFT_FASTA"
SARS_FASTA_PATH = "Sars.txt"
supabase: Client = create_client(SUPABASE_URL, SUPABASE_KEY)
with open(SARS_FASTA_PATH, 'r') as handle:
    wuhan_seq = next(SeqIO.parse(handle, "fasta"))
def align_with_mafft(target_seq: str, wuhan_seq: SeqRecord) -> str:
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f_input:
        records = [wuhan_seq, SeqRecord(Seq(target_seq), id="target", description="")]
        SeqIO.write(records, f_input, "fasta")
        input_path = f_input.name
    with tempfile.NamedTemporaryFile(mode='r', delete=False, suffix='.fasta') as f_output:
        output_path = f_output.name
    subprocess.run(
        ["mafft", "--auto", input_path],
        stdout=open(output_path, 'w'),
        stderr=subprocess.DEVNULL
    )
    alignment = AlignIO.read(output_path, "fasta")
    aligned_seq = str(alignment[1].seq)
    os.remove(input_path)
    os.remove(output_path)
    return aligned_seq

# ---- FETCH RECORDS AND UPDATE ----
def process_and_update(batch_size=500):
    print("Fetching total record count...")
    total_response = supabase.table(TABLE_NAME).select("id", count="exact").execute()
    total_records = total_response.count or 0
    print(f"Total records: {total_records}")

    for start in range(0, total_records, batch_size):
        end = start + batch_size - 1
        print(f"\nFetching records {start} to {end}...")
        response = supabase.table(TABLE_NAME).select("id," + FASTA_COLUMN).range(start, end).execute()
        records = response.data

        for record in records:
            record_id = record["id"]
            print(f"Processing record {record_id}...")
            original_seq = record[FASTA_COLUMN]

            try:
                aligned = align_with_mafft(original_seq, wuhan_seq)
                supabase.table(TABLE_NAME).update({MAFFT_COLUMN: aligned}).eq("id", record_id).execute()
                print(f"✅ Updated record {record_id}")
            except Exception as e:
                print(f"❌ Failed for record {record_id}: {e}")

# ---- RUN ----
if __name__ == "__main__":
    process_and_update()
