# Re-execute the logic since the kernel was reset

from Bio import Entrez, SeqIO
from supabase import create_client
import hashlib

# Setup (replace with your actual keys)
Entrez.email = "wilsonyusda@gmail.com"
SUPABASE_URL = "https://unwqzwvqjdcwmlwafxoo.supabase.co"
SUPABASE_KEY = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InVud3F6d3ZxamRjd21sd2FmeG9vIiwicm9sZSI6ImFub24iLCJpYXQiOjE3NDY4NjU0NzUsImV4cCI6MjA2MjQ0MTQ3NX0.7P3a3cL3JPUidguddJo7HFgy25Ft60UKPpfo4aznqlo"
supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S', 'CCT': 'P', 'CCC': 'P',
    'CCA': 'P', 'CCG': 'P', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
    'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*'
}

def find_all_orfs_fast(dna, codon_table, min_aa_length=30):
    dna = dna.upper()
    stop_codons = {'TAA', 'TAG', 'TGA'}
    orfs = []
    for frame in range(3):
        i = frame
        while i <= len(dna) - 3:
            codon = dna[i:i+3]
            if codon == 'ATG':
                protein = []
                j = i
                while j <= len(dna) - 3:
                    next_codon = dna[j:j+3]
                    if next_codon in stop_codons:
                        break
                    aa = codon_table.get(next_codon)
                    if aa is None:
                        break
                    protein.append(aa)
                    j += 3
                if len(protein) >= min_aa_length:
                    orfs.append(''.join(protein))
                i = j
            else:
                i += 3
    return orfs

def insert_record(row_data):
    try:
        fasta_hash = hashlib.sha256(row_data["FASTA"].encode()).hexdigest()
        row_data["FASTA_HASH"] = fasta_hash

        existing = supabase.table("Nucleotide").select("id").eq("FASTA_HASH", fasta_hash).execute()
        if existing.data:
            print(f"Duplicate: {row_data['accession_id']}")
            return

        supabase.table("Nucleotide").insert(row_data).execute()
        print(f"Inserted: {row_data['accession_id']}")
    except Exception as e:
        print(f"Insert error for {row_data['accession_id']}: {e}")

def process_accessions_in_batches(accession_ids, batch_size=10):
    total_batches = (len(accession_ids) + batch_size - 1) // batch_size
    for batch_num in range(total_batches):
        if batch_num >265:
            start = batch_num * batch_size
            end = start + batch_size
            batch_ids = accession_ids[start:end]
            print(f"Processing batch {batch_num + 1}/{total_batches}: IDs {start + 1}-{min(end, len(accession_ids))}")
            try:
                handle = Entrez.efetch(db="nucleotide", id=",".join(batch_ids), rettype="gb", retmode="text")
                records = list(SeqIO.parse(handle, "genbank"))
                handle.close()
            except Exception as e:
                print(f"Failed fetching batch {batch_num + 1}: {e}")
                continue

            for rec in records:
                location = "N/A"
                collection_date = "N/A"
                for feature in rec.features:
                    if feature.type == "source":
                        location = feature.qualifiers.get('geo_loc_name', ['N/A'])[0]
                        collection_date = feature.qualifiers.get('collection_date', ['N/A'])[0]
                        break

                row_data = {
                    "accession_id": rec.id,
                    "name": rec.name,
                    "length": len(rec.seq),
                    "location": location,
                    "collection_date": collection_date,
                    "FASTA": str(rec.seq),
                    "variants": "Delta",  # forced variant for this dataset
                    "protein_strain": len(find_all_orfs_fast(str(rec.seq), codon_table, min_aa_length=30)),
                }
                insert_record(row_data)

# === Example usage ===
def run_batch_loader():
    combined_ids = set()
    for file_index in range(9, 14):
        file_path = f"{file_index}.acc"
        try:
            with open(file_path, 'r') as f:
                ids = [line.strip() for line in f if line.strip()]
                combined_ids.update(ids)
        except FileNotFoundError:
            print(f"{file_path} not found, skipping.")
    combined_ids = list(combined_ids)
    print(f"Total unique accession IDs: {len(combined_ids)}")
    process_accessions_in_batches(combined_ids, batch_size=10)

run_batch_loader()
