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

    print(f"📥 BAM: {bam}")
    print(f"📄 BED: {bed}")
    print(f"🧬 FASTA: {fasta}")

    # Step 1: Run the main pipeline
    subprocess.run([
        "python", "run/BAM_pipeline_2.py",
        "--bam", bam,
        "--bed", bed,
        "--fasta", fasta
    ], check=True)

    # Step 2: Format outputs for MultiQC
    subprocess.run(["python", "output/BAM_finalize_2.py"], check=True)

    # Step 3: Run MultiQC with your custom config
    subprocess.run([
        "multiqc", ".",
        "-e", "qualimap",
        "-e", "picard",
        "-c", "output/multiqc_config.yaml"
    ], check=True)

    print("✅ QC analysis complete. See 'output/multiqc_report.html' for results.")

if __name__ == "__main__":
    main()
aumoreno@5CD321BNQQ:~/Desktop/BAM_QC$ cat Dockerfile
# Lightweight Python image
FROM python:3.10-slim

# Install system-level dependencies
RUN apt-get update && apt-get install -y \
    samtools \
    openjdk-17-jre \
    git \
    wget \
    && apt-get clean && rm -rf /var/lib/apt/lists/*


# Set working directory
WORKDIR /app

# Copy the entire project into the container
COPY . /app

# Download picard.jar into the run/ directory
RUN wget -O run/picard.jar https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar

# Make scripts executable
RUN chmod +x run/qualimap_v2.3/qualimap \
             run/BAM_pipeline_2.py \
             output/BAM_finalize_2.py \
             run/wrapper.py

# Install Python dependencies (including your custom MultiQC)
RUN pip install --no-cache-dir -r requirements.txt

# Entry point that executes your wrapper
ENTRYPOINT ["python", "run/wrapper.py"]
