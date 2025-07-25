import subprocess
import os
import sys

def main():
    if len(sys.argv) != 7:
        print("Usage: docker run ... bam-qc --bam <file> --bed <file> --fasta <file>")
        sys.exit(1)

    args = dict(zip(sys.argv[1::2], sys.argv[2::2]))
    bam = args.get("--bam")
    bed = args.get("--bed")
    fasta = args.get("--fasta")

    # Validate inputs
    for path, label in zip([bam, bed, fasta], ["BAM", "BED", "FASTA"]):
        if not path or not os.path.isfile(path):
            print(f"Error: {label} file '{path}' not found.")
            sys.exit(1)

    print(f"BAM: {bam}")
    print(f"BED: {bed}")
    print(f"FASTA: {fasta}")

    # Step 1: Run the main pipeline
    subprocess.run([
        "python", "run/BAM_pipeline_2.py",
        "--bam", bam,
        "--bed", bed,
        "--fasta", fasta
    ], check=True)

    # Step 2: Format outputs for MultiQC
    subprocess.run(["python", "output/BAM_finalize_2.py"], check=True)


if __name__ == "__main__":
    main()
