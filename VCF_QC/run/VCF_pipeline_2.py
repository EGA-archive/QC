#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Process a VCF or VCF.GZ file to extract QC metrics and summary info.")
    parser.add_argument("--vcf", required=True, help="Path to the input VCF or VCF.GZ file")
    args = parser.parse_args()

    vcf_path = os.path.abspath(args.vcf)
    if not os.path.exists(vcf_path):
        print(f"[ERROR] File not found: {vcf_path}")
        sys.exit(1)

    if vcf_path.endswith(".vcf.gz"):
        base = os.path.basename(vcf_path).replace(".vcf.gz", "")
        vcf_cat = f"zcat {vcf_path}"
    elif vcf_path.endswith(".vcf"):
        base = os.path.basename(vcf_path).replace(".vcf", "")
        vcf_cat = f"cat {vcf_path}"
    elif vcf_path.endswith(".bcf"):
        base = os.path.basename(vcf_path).replace(".vcf", "")
        vcf_cat = f"cat {vcf_path}"
    elif vcf_path.endswith(".bcf.gz"):
        base = os.path.basename(vcf_path).replace(".vcf", "")
        vcf_cat = f"cat {vcf_path}"
    else:
        print("[ERROR] Input file must end in .vcf, .bcf, .vcf.gz or bcf.gz")
        sys.exit(1)

    # Output folder: one level above script, in ../output
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.abspath(os.path.join(script_dir, "..", "output"))
    os.makedirs(output_dir, exist_ok=True)

    print(f"[INFO] Processing VCF file: {vcf_path}")
    print(f"[INFO] Saving outputs to: {output_dir}")

    # Output files
    index_file = os.path.join(output_dir, f"{base}.csi")
    random_snps_file = os.path.join(output_dir, f"random_snps.txt")
    sample_list_file = os.path.join(output_dir, f"samples.txt")
    chrom_list_file = os.path.join(output_dir, f"chromosomes.txt")
    vcf_version_file = os.path.join(output_dir, f"version.txt")
    header_file = os.path.join(output_dir, f"header.txt")
    stats_sample_file = os.path.join(output_dir, f"stats_sample.txt")
    stats_file = os.path.join(output_dir, f"stats.txt")
    # Step 1: Index the VCF
    print("[1] Indexing the VCF...")
    try:
        subprocess.run(f"bcftools index -o {index_file} {vcf_path}", shell=True, check=True)
    except subprocess.CalledProcessError:
        print("[ERROR] Failed to index the VCF.")

    # Step 2: Extract 1000 random rsID SNPs using zcat or cat
    print("[2] Extracting 1000 random rsID SNPs...")
    try:
        cmd = f"{vcf_cat} | cut -f1,2,3 | awk '$3 ~ /rs/' | shuf -n 1000 > {random_snps_file}"
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        print("[ERROR] Failed to extract SNPs.")


    # Step 3: Sample list
    print("[3] Extracting sample list...")
    try:
        subprocess.run(f"bcftools query -l {vcf_path} > {sample_list_file}", shell=True, check=True)
    except subprocess.CalledProcessError:
        print("[ERROR] Failed to extract sample list.")


    # Step 4: Unique chromosomes
    print("[4] Extracting unique chromosomes...")
    try:
        subprocess.run(f"bcftools query -f '%CHROM\\n' {vcf_path} | uniq > {chrom_list_file}", shell=True, check=True)
    except subprocess.CalledProcessError:
        print("[ERROR] Failed to extract chromosomes.")


    # Step 5: VCF version (from first line)
    print("[5] Detecting VCF version...")
    try:
        cmd = f"{vcf_cat} | head -n 1"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        with open(vcf_version_file, "w") as out:
            first_line = result.stdout.strip()
            if "VCF" in first_line:
                version = first_line.split("VCF")[1].strip()
                out.write(version + "\n")
                print(f"    → Detected version: {version}")
            else:
                out.write("Version not found\n")
                print("    → No version string found.")
    except subprocess.CalledProcessError:
        print("[ERROR] Could not read VCF version.")


    # Step 6: Extract header excluding #CHROM
    print("[6] Extracting header (excluding #CHROM)...")
    try:
        cmd = f"{vcf_cat} | grep '^#' | grep -vi '^#CHROM' > {header_file}"
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        print("[ERROR] Failed to extract header.")



    # Step 7: Per sample statistics
    print("[7] Calculating bcftools sample stats...")
    try:
        cmd = f"bcftools stats -S {sample_list_file} {vcf_path} > {stats_sample_file}"
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        print("[ERROR] Failed to generate stats.")


    # Step 8: Statistics
    print("[8] Calculating bcftools stats...")
    try:
        cmd = f"bcftools stats {vcf_path} > {stats_file}"
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        print("[ERROR] Failed to generate stats.")


    # Final summary
    print("\n✅ VCF processing completed.")
    print("Generated files:")
    for f in [
        index_file, random_snps_file, sample_list_file, chrom_list_file,
        vcf_version_file, header_file, stats_file, stats_sample_file
    ]:
        print(f" - {f}")

if __name__ == "__main__":
    main()
