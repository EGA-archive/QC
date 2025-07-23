import subprocess
import os
import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: docker run ... bam-qc /path/to/input.bam")
        sys.exit(1)

    bam_path = sys.argv[1]
    if not os.path.isfile(bam_path):
        print(f"Error: BAM file '{bam_path}' not found inside container.")
        sys.exit(1)

    print(f"📥 Received BAM file: {bam_path}")

    # Step 1: Run the main pipeline
    subprocess.run(["python", "run/BAM_pipeline_2.py", bam_path], check=True)

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
