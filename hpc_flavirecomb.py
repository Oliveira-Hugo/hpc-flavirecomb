#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

def run(cmd, cwd=None, log=None):
    """Run a command, print it, and optionally log output."""
    print(f"\n=== Running: {cmd} ===\n")
    result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)

    if log:
        log.write(f"\n--- COMMAND: {cmd} ---\n")
        log.write(result.stdout)
        log.write(result.stderr)

    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(f"Command failed: {cmd}")

    return result.stdout


def main():
    project_root = Path(__file__).resolve().parent

    TREE_PIPE = project_root / "tree_pipeline" / "scripts" / "tree_pipeline.py"
    ALT_PIPE_DIR = project_root / "alternative_trees_pipeline" / "scripts"
    SNP_PIPE_DIR = project_root / "snp_pipeline" / "scripts"

    input_fasta = project_root / "input" / "annotated_denv_genomes_nm.fasta"

    with open(project_root / "workflow.log", "w") as LOG:

        LOG.write("=== Starting Workflow ===\n")

        # 1) TREE PIPELINE (IQ-TREE → TNT input)
        LOG.write("\n\n### TREE PIPELINE ###\n")
        
        alignment_file = project_root / "tree_pipeline" / "data" / "alignment.fasta"
        tree_outdir = project_root / "tree_pipeline" / "results"

        if TREE_PIPE.exists():
            if not alignment_file.exists():
                raise FileNotFoundError(f"Tree alignment missing: {alignment_file}")
        
            run(
                f"python3 {TREE_PIPE} "
                f"--alignment {alignment_file} "
                f"--outdir {tree_outdir} "
                f"--prefix mytree "
                f"--threads 8",
                log=LOG
            )
        else:
            LOG.write("Tree pipeline script NOT FOUND — skipping\n")


        # 2) ALTERNATIVE TREES PIPELINE
        LOG.write("\n\n### ALTERNATIVE TREES PIPELINE ###\n")

        # 2.1 generate alt alignments
        alt00 = ALT_PIPE_DIR / "00_prepare_alt_alignments.py"
        if alt00.exists():
            run(f"python3 {alt00}", log=LOG)

        # 2.2 generate TNT scripts
        alt01 = ALT_PIPE_DIR / "01_prepare_tnt_scripts.py"
        if alt01.exists():
            run(f"python3 {alt01}", log=LOG)

        # 2.3 run TNT batch jobs
        alt02 = ALT_PIPE_DIR / "02_run_tnt_scripts.sh"
        if alt02.exists():
            alt02.chmod(0o755)
            run(str(alt02), cwd=ALT_PIPE_DIR, log=LOG)

        # 3) SNP PIPELINE
        LOG.write("\n\n### SNP PIPELINE ###\n")

        # 3.1 Run Snipit automatically
        snp00 = SNP_PIPE_DIR / "00_auto_snipit.py"
        if snp00.exists():
            run(f"python3 {snp00}", log=LOG)

        # 3.2 Snipit splitting (bash script)
        snp01 = SNP_PIPE_DIR / "01_run_snipit.sh"
        if snp01.exists():
            snp01.chmod(0o755)
            run(str(snp01), cwd=SNP_PIPE_DIR, log=LOG)

        # 3.3 Calculate SNPs
        snp02 = SNP_PIPE_DIR / "02_calculate_snps.py"
        if snp02.exists():
            run(f"python3 {snp02}", log=LOG)


        LOG.write("\n=== WORKFLOW COMPLETED SUCCESSFULLY ===\n")

    print("\nSee workflow.log for details.\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("\n Execution failed:", e)
        sys.exit(1)

