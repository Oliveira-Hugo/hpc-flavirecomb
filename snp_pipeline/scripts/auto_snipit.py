#!/usr/bin/env python3
from Bio import SeqIO
import csv
import os
from pathlib import Path

root = Path(__file__).resolve().parent.parent

fasta_file = root / "data" / "annotated_denv_genomes_nm.fasta"
table_file = root / "data" / "recomb_and_parents.csv"
output_dir = root / "data" / "fragments"
alignment_length = 15213  # fixed alignment length

# ==== CREATE OUTPUT DIR ====
os.makedirs(output_dir, exist_ok=True)

# ==== LOAD FASTA ====
seq_dict = {rec.id: rec for rec in SeqIO.parse(fasta_file, "fasta")}

# ==== READ TABLE AND PROCESS ====
with open(table_file, newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=",")
    for row in reader:
        recombinant = row["Recombinant"].strip()
        minor = row["Minor parent"].strip()
        major = row["Major parent"].strip()
        begin = int(row["Begin"])
        end = int(row["End"])

        # Filter out "Unknown"
        taxa = [t for t in [recombinant, minor, major] if t != "Unknown"]

        # Skip if no taxa found
        if not taxa:
            continue

        # Check sequences exist
        missing = [t for t in taxa if t not in seq_dict]
        if missing:
            print(f"Warning: Missing sequences for {missing}, skipping line.")
            continue

        # Output base name
        base_name = f"snipit_{recombinant}_{minor}_{major}".replace("Unknown", "NA")
        out_fasta_path = os.path.join(output_dir, f"{base_name}.fasta")

        # Write original sequences
        SeqIO.write([seq_dict[t] for t in taxa], out_fasta_path, "fasta")

        # Create Fragment 1 (Begin to End)
        frag1_records = []
        for t in taxa:
            seq = seq_dict[t].seq
            frag1_seq = seq[begin - 1:end]  # Python is 0-based
            frag1_records.append(seq_dict[t][:0])  # clone empty
            frag1_records[-1].seq = frag1_seq

        frag1_path = os.path.join(output_dir, f"{base_name}_frag1_{begin}_{end}.fasta")
        SeqIO.write(frag1_records, frag1_path, "fasta")

        # Create Fragment 2 (1 to Begin-1 + End to alignment_length)
        frag2_records = []
        for t in taxa:
            seq = seq_dict[t].seq
            frag2_seq = seq[0:begin - 1] + seq[end:alignment_length]
            frag2_records.append(seq_dict[t][:0])  # clone empty
            frag2_records[-1].seq = frag2_seq

        frag2_path = os.path.join(output_dir, f"{base_name}_frag2_outside_{begin}_{end}.fasta")
        SeqIO.write(frag2_records, frag2_path, "fasta")

print("Done.")
