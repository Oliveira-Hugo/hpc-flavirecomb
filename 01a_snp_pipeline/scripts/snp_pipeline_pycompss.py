from pathlib import Path
from pycompss.api.task import task
from pycompss.api.binary import binary
from pycompss.api.parameter import FILE_IN
from pycompss.api.api import compss_barrier
import subprocess
import sys

@binary(
    binary="/scratch/cenapadrjsd/hugo.oliveira2/hpc-flavirecomb/01a_snp_pipeline/scripts/01_run_snipit.sh"
)
@task(fasta=FILE_IN)
def run_snipit(fasta, out_dir):
    pass

def generate_fragments():
    """Executa o script 00_auto_snipit.py no master para criar os FASTAs."""
    script_path = Path("/scratch/cenapadrjsd/hugo.oliveira2/hpc-flavirecomb/01a_snp_pipeline/scripts/00_auto_snipit.py")
    result = subprocess.run([sys.executable, str(script_path)], check=True)
    if result.returncode != 0:
        raise RuntimeError("Falha ao gerar fragmentos FASTA com 00_auto_snipit.py")

def calculate_snps():
    """Executa o script 02_calculate_snps.py no master para gerar estatísticas de SNPs."""
    script_path = Path("/scratch/cenapadrjsd/hugo.oliveira2/hpc-flavirecomb/01a_snp_pipeline/scripts/02_calculate_snps.py")
    result = subprocess.run([sys.executable, str(script_path)], check=True)
    if result.returncode != 0:
        raise RuntimeError("Falha ao calcular SNPs com 02_calculate_snps.py")

def main():
    print("MASTER: starting")
    
    print("MASTER: generating fragment FASTAs")
    generate_fragments()

    frag_dir = Path("/scratch/cenapadrjsd/hugo.oliveira2/hpc-flavirecomb/01a_snp_pipeline/results/fragments")
    out_root = Path("/scratch/cenapadrjsd/hugo.oliveira2/hpc-flavirecomb/01a_snp_pipeline/results/snipit_outputs").resolve()
    out_root.mkdir(exist_ok=True)

    fasta_files = sorted(frag_dir.glob("*.fasta"))
    print(f"MASTER: found {len(fasta_files)} FASTA files")

    for fasta in fasta_files:
        task_dir = out_root / fasta.stem
        task_dir.mkdir(exist_ok=True)

        print(f"MASTER: scheduling {fasta.name}")
        run_snipit(str(fasta), str(task_dir))

    print("MASTER: waiting for tasks")
    compss_barrier()
    print("MASTER: Snipit tasks completed")

    print("MASTER: calculating SNP statistics")
    calculate_snps()
    print("MASTER: done")

if __name__ == "__main__":
    main()
