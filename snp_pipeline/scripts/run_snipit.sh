#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda activate snipit-env

FRAG_DIR="../data/fragments"
cd "$FRAG_DIR" || { echo "Cannot access $FRAG_DIR"; exit 1; }

for fasta_file in *.fasta; do
    if [[ -f "$fasta_file" ]]; then
        echo "Processing $fasta_file..."

        snipit "$fasta_file" -s

        if [[ -f "snps.csv" ]]; then
            new_name="${fasta_file%.fasta}.csv"
            mv snps.csv "$new_name"
            echo "Renamed snps.csv to $new_name"
        else
            echo "Warning: snps.csv not found for $fasta_file"
        fi
    fi
done

