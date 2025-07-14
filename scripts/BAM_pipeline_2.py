#!/usr/bin/env python3
import subprocess
import argparse
import os
import pysam
import json
import math
import statistics
from collections import defaultdict, Counter

def run_command(command):
    print(f"[INFO] Running command: {command}")
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error executing the line: {command}")
        print(result.stderr)
    else:
        print(result.stdout)

def convert_cram_to_bam(cram_path):
    bam_path = cram_path.replace(".cram", ".bam")
    print(f"Converting {cram_path} to BAM: {bam_path}")
    cmd = f"samtools view -b -o {bam_path} {cram_path}"
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise RuntimeError(f"Error converting {cram_path} to BAM.")
    return bam_path

def create_bam_index(bam_path):
    bai_path = bam_path + ".bai"
    if not os.path.exists(bai_path):
        print(f" {os.path.basename(bai_path)} was not found. Creating index with pysam.index()…")
        subprocess.run(["samtools", "index", bam_path], check=True)


def parse_bam_header(bam_file, output_file):
    result = subprocess.run(["samtools", "view", "-H", bam_file], capture_output=True, text=True)
    header_lines = result.stdout.split('\n')
    data = {
        "Version": "N/A",
        "Sort Order": "N/A",
        "Platform": "N/A",
        "Tools": set(),
        "Genome Assembly Identifier": "N/A",
        "Sample" : "Sample"
    }
    for line in header_lines:
        fields = line.strip().split('\t')
        if line.startswith("@HD"):
            for field in fields:
                if field.startswith("VN:"):
                    data["Version"] = field.split(":")[1]
                elif field.startswith("SO:"):
                    data["Sort Order"] = field.split(":")[1]
        elif line.startswith("@RG"):
            for field in fields:
                if field.startswith("PL:"):
                    data["Platform"] = field.split(":")[1]
                elif field.startswith("SM:"):
                    data["Sample"] = field.split(":")[1]
        elif line.startswith("@PG"):
            for field in fields:
                if field.startswith("ID:") or field.startswith("PN:"):
                    data["Tools"].add(field.split(":")[1])
        elif line.startswith("@SQ"):
            for field in fields:
                if field.startswith("AS:"):
                    data["Genome Assembly Identifier"] = field.split(":")[1]
    data["Tools"] = sorted(data["Tools"])
    with open(output_file, "w") as f:
        f.write(f"Version\t{data['Version']}\n")
        f.write(f"Sort Order\t{data['Sort Order']}\n")
        f.write(f"Platform\t{data['Platform']}\n")
        f.write(f"Genome Assembly Identifier\t{data['Genome Assembly Identifier']}\n")
        f.write(f"Sample\t{data['Sample']}\n")
        f.write(f"Tools\t{' | '.join(data['Tools'])}\n")
    print(f"Header info saved to {output_file}")

def calculate_stat_from_counter(counter):
    values = []
    for cov, freq in counter.items():
        values.extend([cov] * freq)
    if not values:
        return 0.0
    return round(statistics.median(values), 2)

def calculate_coverage_stat_streaming(bam_file, window_size, output_file, alias_json_path="chrnames.json"):
    with open(alias_json_path, "r") as f:
        alias_map = json.load(f)

    create_bam_index(bam_file)
    bam = pysam.AlignmentFile(bam_file, "rb")

    ref_id_to_name = {}
    cumulative_offsets = {}
    ref_name_to_length = {}
    offset = 0
    matched_any = False

    for i in range(bam.nreferences):
        ref = bam.get_reference_name(i)
        canonical = alias_map.get(ref)
        if canonical:
            matched_any = True
            length = bam.get_reference_length(ref)
            ref_id_to_name[i] = canonical
            ref_name_to_length[canonical] = length
            if canonical not in cumulative_offsets:
                cumulative_offsets[canonical] = offset
                offset += length

    if not matched_any:
        print("No valid names in the alias. Using fallback.")
        canonical_names = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
        for i, canonical in enumerate(canonical_names):
            if i >= bam.nreferences:
                break
            ref = bam.references[i]
            length = bam.get_reference_length(ref)
            ref_id_to_name[i] = canonical
            ref_name_to_length[canonical] = length
            cumulative_offsets[canonical] = offset
            offset += length
        ref_map = dict(zip(bam.references, canonical_names))
    else:
        ref_map = {ref: alias_map[ref] for ref in bam.references if ref in alias_map}

    coverage_by_window = defaultdict(Counter)
    chrM_coverages = []

    proc = subprocess.Popen(["samtools", "depth", "-a", bam_file], stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        ref, pos, cov = line.strip().split("\t")
        if ref not in ref_map:
            continue
        canonical = ref_map[ref]
        pos = int(pos)
        cov = int(cov)
        if canonical == "chrM":
            chrM_coverages.append(cov)
            continue
        global_pos = cumulative_offsets[canonical] + pos
        window_start = (global_pos // window_size) * window_size
        coverage_by_window[window_start][cov] += 1
    proc.stdout.close()
    proc.wait()

    offset_to_chrom = sorted(cumulative_offsets.items(), key=lambda x: x[1])
    def get_chrom(global_pos):
        for i in range(len(offset_to_chrom) - 1):
            chrom, start = offset_to_chrom[i]
            _, next_start = offset_to_chrom[i + 1]
            if start <= global_pos < next_start:
                return chrom
        return offset_to_chrom[-1][0]

    total_length = sum(ref_name_to_length[c] for c in ref_name_to_length if c != "chrM")
    num_windows = math.ceil(total_length / window_size)
    result = []

    for i in range(num_windows):
        global_start = i * window_size
        chrom = get_chrom(global_start)
        relative_start = global_start
        counter = coverage_by_window.get(global_start, Counter())
        cov_value = calculate_stat_from_counter(counter, mode) if counter else 0.0
        result.append([chrom, relative_start, cov_value])

    if chrM_coverages:
        artificial_pos = num_windows * window_size
        chrM_stat = round(statistics.median(chrM_coverages), 2)
        result.append(["chrM", artificial_pos, chrM_stat])
    else:
        print(" [WARNING] chrM does not have a coverage")

    with open(output_file, "w") as f:
        json.dump(result, f, separators=(",", ":"))

    print(f"[INFO] Result saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Análisis de archivos BAM o CRAM con herramientas bioinformáticas.")
    parser.add_argument("--bam", required=True, help="Archivo BAM o CRAM de entrada")
    parser.add_argument("--bed", required=True, help="Archivo BED necesario para algunas herramientas")
    parser.add_argument("--mode", choices=["mean", "median"], default="mean", help="Modo de cálculo: media o mediana")
    parser.add_argument("--fasta", required= True, help="Archivo FASTA")
    args = parser.parse_args()

    input_path = os.path.abspath(args.bam)
    bed_file = os.path.abspath(args.bed)
    fasta_file=os.path.abspath(args.fasta)

    # If the file is a .cram, convert it to .bam
    if input_path.endswith(".cram"):
        bam_file = convert_cram_to_bam(input_path)
    elif input_path.endswith(".bam"):
        bam_file = input_path
    else:
        raise ValueError("Input file has to be BAM or CRAM")

    output_dir = os.path.join(os.path.dirname(bam_file), "bam_analysis_results")
    os.makedirs(output_dir, exist_ok=True)
    alias_json = os.path.join(os.path.dirname(__file__), "chrnames.json")

    header_output_file = os.path.join(output_dir, "bam_header_info.txt")
    parse_bam_header(bam_file, header_output_file)

    run_command(f"/bio-scratch/angel/qualimap-install-attempt/qualimap_v2.3/qualimapACT bamqc -bam {bam_file} -c --outdir {output_dir} --java-mem-size=16G")
    run_command(f"read_distribution.py -i {bam_file} -r {bed_file} > {output_dir}/read_distribution.txt")
    run_command(f"java -jar /bio-scratch/Aurora/BAM_report/definitivo/picard.jar  CollectAlignmentSummaryMetrics           R= {fasta_file}           I= {bam_file}          O={output_dir}/picard_output.txt")
    run_command(f"java -jar /bio-scratch/Aurora/BAM_report/definitivo/picard.jar CollectQualityYieldMetrics I= {bam_file} O= {output_dir}/collect_bases_metrics.txt")

    window_cov = 3_000_000
    cov_output_file = os.path.join(output_dir, f"coverage_median_{window_cov // 1_000_000}Mb_complete.json")
    calculate_coverage_stat_streaming(bam_file, window_cov, cov_output_file, alias_json_path=alias_json)

if __name__ == "__main__":
    main()
