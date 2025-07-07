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
    print(f"Ejecutando: {command}")
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error ejecutando el comando: {command}")
        print(result.stderr)
    else:
        print(result.stdout)

def convert_cram_to_bam(cram_path):
    bam_path = cram_path.replace(".cram", ".converted.bam")
    print(f"🌀 Convirtiendo {cram_path} a BAM: {bam_path}")
    cmd = f"samtools view -b -o {bam_path} {cram_path}"
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise RuntimeError(f"❌ Error al convertir {cram_path} a BAM.")
    return bam_path

def create_bam_index(bam_path):
    bai_path = bam_path + ".bai"
    if not os.path.exists(bai_path):
        print(f"🔧 No se encontró {os.path.basename(bai_path)}. Creando índice con pysam.index()…")
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
def median_from_counter(counter):
    total = sum(counter.values())
    if total == 0:
        return 0.0
    mid = total // 2
    acc = 0
    sorted_items = sorted(counter.items())
    for i, (value, count) in enumerate(sorted_items):
        acc += count
        if acc > mid:
            return float(value)
        elif acc == mid and total % 2 == 0:
            for j in range(i + 1, len(sorted_items)):
                return (value + sorted_items[j][0]) / 2.0
            return float(value)
    return 0.0
def calculate_mapq_median(bam_file, window_size, output_file, alias_json_path="chrnames.json"):
    with open(alias_json_path, "r") as f:
        alias_map = json.load(f)
    create_bam_index(bam_file)
    bam = pysam.AlignmentFile(bam_file, "rb")
    mapq_by_global_window = defaultdict(list)
    ref_id_to_name = {}
    cumulative_offsets = {}
    ref_name_to_length = {}
    offset = 0
    matched_any = False
    for i in range(bam.nreferences):
        original_name = bam.get_reference_name(i)
        canonical = alias_map.get(original_name)
        if canonical:
            matched_any = True
            length = bam.get_reference_length(original_name)
            ref_id_to_name[i] = canonical
            ref_name_to_length[canonical] = length
            if canonical not in cumulative_offsets:
                cumulative_offsets[canonical] = offset
                offset += length
    if not matched_any:
        print("No se encontraron nombres válidos en el alias. Usando fallback forzado.")
        canonical_names = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
        for i, canonical in enumerate(canonical_names):
            if i >= bam.nreferences:
                print(f"Faltan contigs para cubrir hasta {canonical}")
                break
            length = bam.get_reference_length(i)
            ref_id_to_name[i] = canonical
            ref_name_to_length[canonical] = length
            cumulative_offsets[canonical] = offset
            offset += length
    valid_ref_ids = set(ref_id_to_name.keys())
    for read in bam:
        if read.reference_id not in valid_ref_ids or read.reference_start < 0:
            continue
        chrom = ref_id_to_name[read.reference_id]
        global_pos = cumulative_offsets[chrom] + read.reference_start
        window_start = (global_pos // window_size) * window_size
        mapq_by_global_window[window_start].append(read.mapping_quality)
    bam.close()
    total_length = sum(ref_name_to_length.values())
    num_windows = math.ceil(total_length / window_size)
    offset_to_chrom = sorted(cumulative_offsets.items(), key=lambda x: x[1])
    def get_chrom_for_pos(global_pos):
        for i in range(len(offset_to_chrom) - 1):
            chrom, start = offset_to_chrom[i]
            next_chrom, next_start = offset_to_chrom[i + 1]
            if start <= global_pos < next_start:
                return chrom
        return offset_to_chrom[-1][0]
    result = []
    for i in range(num_windows):
        global_start = i * window_size
        chrom = get_chrom_for_pos(global_start)
        values = mapq_by_global_window.get(global_start, [])
        median_mapq = round(statistics.median(values), 2) if values else 0.0
        result.append([chrom, global_start, median_mapq])
    chroms_included = {entry[0] for entry in result}
    chrM_aliases = [k for k, v in alias_map.items() if v == "chrM"]
    mt_offset = None
    if "chrM" not in chroms_included:
        for alias in chrM_aliases:
            if alias in cumulative_offsets:
                mt_offset = cumulative_offsets[alias]
                break
    if mt_offset is not None:
        mapqs = []
        with pysam.AlignmentFile(bam_file, "rb") as bam_mt:
            for alias in chrM_aliases:
                if alias in bam_mt.references:
                    for r in bam_mt.fetch(alias):
                        if not r.is_unmapped:
                            mapqs.append(r.mapping_quality)
                    break
        if mapqs:
            median_mapq = round(statistics.median(mapqs), 2)
            result.append(["chrM", mt_offset, median_mapq])
            print(f":white_check_mark: chrM añadido con mediana MAPQ={median_mapq}")
        else:
            print(":warning: chrM sin lecturas mapeadas")
    with open(output_file, "w") as f:
        json.dump(result, f, separators=(",", ":"))
    print(f"Mediana de MAPQ calculada y guardada en {output_file}")


def calculate_stat_from_counter(counter, mode="mean"):
    values = []
    for cov, freq in counter.items():
        values.extend([cov] * freq)
    if not values:
        return 0.0
    return round(statistics.mean(values), 2) if mode == "mean" else round(statistics.median(values), 2)

def calculate_coverage_stat_streaming(bam_file, window_size, output_file, alias_json_path="chrnames.json", mode="mean"):
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
        print("⚠️ No se encontraron nombres válidos en el alias. Usando fallback.")
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
        chrM_stat = round(statistics.mean(chrM_coverages), 2) if mode == "mean" else round(statistics.median(chrM_coverages), 2)
        result.append(["chrM", artificial_pos, chrM_stat])
        print(f":white_check_mark: chrM añadido con {mode} = {chrM_stat}")
    else:
        print(":warning: chrM no tiene cobertura")

    with open(output_file, "w") as f:
        json.dump(result, f, separators=(",", ":"))

    print(f"✅ Resultado guardado en {output_file}")

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

    # Si el archivo es .cram, conviértelo a .bam
    if input_path.endswith(".cram"):
        bam_file = convert_cram_to_bam(input_path)
    elif input_path.endswith(".bam"):
        bam_file = input_path
    else:
        raise ValueError("❌ El archivo de entrada debe ser .bam o .cram")

    output_dir = os.path.join(os.path.dirname(bam_file), "bam_analysis_results")
    os.makedirs(output_dir, exist_ok=True)
    alias_json = os.path.join(os.path.dirname(__file__), "chrnames.json")

    header_output_file = os.path.join(output_dir, "bam_header_info.txt")
    parse_bam_header(bam_file, header_output_file)

    run_command(f"/bio-scratch/angel/qualimap-install-attempt/qualimap_v2.3/qualimapACT bamqc -bam {bam_file} -c --outdir {output_dir} --java-mem-size=16G")
    run_command(f"read_distribution.py -i {bam_file} -r {bed_file} > {output_dir}/read_distribution.txt")
    run_command(f"java -jar /bio-scratch/Aurora/BAM_report/definitivo/picard.jar  CollectAlignmentSummaryMetrics           R= {fasta_file}           I= {bam_file}          O=picard_output.txt")
    run_command(f"java -jar /bio-scratch/Aurora/BAM_report/definitivo/picard.jar QualityYieldMetrics I= {bam_file} O= collect_bases_metrics.txt")
    run_command(f"samtools view -H {bam_file} > header.txt")

    window_mapq = 3_000_000
    mapq_output_file = os.path.join(output_dir, f"mapq_median_{window_mapq // 1_000_000}Mb_complete.json")
    calculate_mapq_median(bam_file, window_mapq, mapq_output_file, alias_json_path=alias_json)

    window_cov = 3_000_000
    cov_output_file = os.path.join(output_dir, f"coverage_median_{window_cov // 1_000_000}Mb_complete.json")
    calculate_coverage_stat_streaming(bam_file, window_cov, cov_output_file, alias_json_path=alias_json)

if __name__ == "__main__":
    main()
