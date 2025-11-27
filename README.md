# hpc-flavirecomb

This project develops a scalable bioinformatics workflow designed to detect, analyze, and characterize recombination events in flavivirus genomes. The workflow designed to be compatible with High Performance Computing (HPC) environments. 

hpc-flavirecomb integrates sequence alignemnt and annotation, recombination detection, SNP profiling, phylogenetic inference, and evaluation of recombination impact on phylogenetic tree.

## Steps

01_hpc_flavirecomb/  
├── 00_FLAVi/         # Annotation, quality control, and MSA generation for flavivirus genomes  
├── 01_RDP5/          # Recombination event detection (GUI-only, to be replaced by OpenRDP or related tool)  
├── 02_SNIPIT/        # SNP identification on recombinant and parental sequences  
├── 03_IQTREE_TNT/    # Phylogenetic tree inference using IQ-TREE2 and TNT    
└── 04_TNT_YBYRA/         # Statistical summary of SNPs, branch lengths, phylogenetic distances between recombinants, and recombinants impact on tree topology

## Current Limitations

Step 00 (FLAVi) is not yet implemented because it depends on the output structure of Step 01.

Step 01 (RDP5) is currently incompatible with HPC environments because RDP5 is GUI-only.
A replacement using OpenRDP (CLI-capable) is under evaluation and will be integrated in future versions.

Step 04 still not implement YBYRÁ for comparing alternative trees.

## Dependencies:
1. Python (≥3.9)  
-biopython  
-pandas  
-snipit  

2. IQ-TREE v2.0.7

3. TNT v1.6

4. YBYRÁ
-svgwrite
