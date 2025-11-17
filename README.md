# hpc-flavirecomb

This project develops a scalable bioinformatics workflow designed to detect, analyze, and characterize recombination events in flavivirus genomes. The workflow designed to be compatible with High Performance Computing (HPC) environments. 

hpc-flavirecomb integrates sequence alignemnt and annotation, recombination detection, SNP profiling, phylogenetic inference, and evaluation of recombination impact on phylogenetic tree.

## Steps

01_hpc_flavirecomb/  
├── 00_FLAVi/         # Annotation, quality control, and MSA generation for flavivirus genomes  
├── 01_RDP5/          # Recombination event detection (GUI-only, to be replaced by OpenRDP)  
├── 02_SNIPIT/        # SNP extraction from recombinant and parental sequences  
├── 03_IQTREE_TNT/    # Phylogenetic tree inference using IQ-TREE3 and TNT  
├── 04_YBYRA/         # Quantification of recombinant impact on tree topology/branch lengths  
└── 05_STATISTICS/    # Statistical summary of SNPs, branch lengths, and phylogenetic distances

## Current Limitations

Step 00 (FLAVi) is not yet implemented because it depends on the output structure of Step 01.

Step 01 (RDP5) is currently incompatible with HPC environments because RDP5 is GUI-only.
A replacement using OpenRDP (CLI-capable) is under evaluation and will be integrated in future versions.

## Dependencies:
1. Python (≥3.9)
-biopython
-pandas
-snipit

2. IQ-TREE 3

3. TNT v1.6
