from Bio import SeqIO
import os

input_fasta = "/home/hugo/hpc_flavirecomb/01b_tree_pipeline/data/alignment.fasta"
alternative_alignments_dir = "/home/hugo/hpc_flavirecomb/01c_alternative_trees_pipeline/results"
tnt_scripts_dir = alternative_alignments_dir
os.makedirs(tnt_scripts_dir, exist_ok=True)

sequence_ids = [seq.id for seq in SeqIO.parse(input_fasta, "fasta")]

tnt_template = """log tnt_{terminal}.log ;
sect : slack 10 ;
mxram 2000 ;
nstates dna ;
taxname +100 ;
proc {alignment_file} ;
hold 10000 ;
xmult = level 3 chklevel 5 hits 100 rep 1000 ;
best ;
length ;
tsave mpts_{terminal}.tnt ;
taxname = ;
tplot ;
save ;
tsave / ;
tsave * consensus_{terminal}.tnt ;
nelsen * ;
save / ;
tsave / ;
CLS ;
quit ;
"""

count = 0
for terminal in sequence_ids:
    terminal_safe = terminal.replace(".", "_").replace("-", "_")

    alignment_file_nexus = os.path.join(
        alternative_alignments_dir,
        f"alternative_alignment_{terminal_safe}_removed.nexus"
    )

    if not os.path.isfile(alignment_file_nexus):
        print(f"Warning: Nexus file not found for {terminal}: {alignment_file_nexus}")
        continue

    alignment_file_nexus = os.path.abspath(alignment_file_nexus).strip()

    script_content = tnt_template.format(
        terminal=terminal_safe,
        alignment_file=alignment_file_nexus
    )

    script_path = os.path.join(tnt_scripts_dir, f"script_{terminal_safe}.RUN")
    with open(script_path, "w", newline="\n") as f:
        f.write(script_content)

    count += 1

print(f"Generated {count} TNT scripts in {tnt_scripts_dir}")
