from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import os

input_fasta = "/home/hugo/hpc_flavirecomb/tree_pipeline/data/alignment.fasta"
output_dir = "/home/hugo/hpc_flavirecomb/alternative_trees_pipeline/results"
os.makedirs(output_dir, exist_ok=True)

sequences = list(SeqIO.parse(input_fasta, "fasta"))

for seq in sequences:
    seq_id_safe = seq.id.replace(".", "_")
    
    output_nexus = os.path.join(output_dir, f"alternative_alignment_{seq_id_safe}_removed.nexus")
    
    filtered_seqs = [s for s in sequences if s.id != seq.id]
    
    dna_seqs = [
        SeqRecord(
            Seq(str(s.seq)),
            id=s.id.replace(".", "_"),           # safe ID
            description="",
            annotations={"molecule_type": "DNA"}
        )
        for s in filtered_seqs
    ]
    
    alignment = MultipleSeqAlignment(dna_seqs, annotations={"molecule_type": "DNA"})
    
    SeqIO.write(alignment, output_nexus, "nexus")

print(f"Generated {len(sequences)} alternative NEXUS alignments in {output_dir}")
