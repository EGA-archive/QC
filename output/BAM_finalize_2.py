#!/usr/bin/env python3

import os
import re
import pandas as pd
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

# === PART 1: Process genome_results.txt, picard_output.txt, collect_bases_metrics.txt  ===

def parse_genome_results(file_path):
    patterns = {
        "number of reads": (r"number of reads\s*=\s*([\d,]+)", "Number of reads"),
        "number of mapped reads": (r"number of mapped reads\s*=\s*([\d,]+)", "Number of mapped reads"),
        "number of mapped paired reads \(both in pair\)": (r"number of mapped paired reads \(both in pair\)\s*=\s*([\d,]+)", "Number of mapped paired reads (both in pair)"),
        "number of mapped paired reads \(singletons\)": (r"number of mapped paired reads \(singletons\)\s*=\s*([\d,]+)", "Number of mapped paired reads (singletons)"),
        "Median insert size": (r"median insert size\s*=\s*([\d,.]+)", "Median insert size"),
        "Mean mapping quality": (r"mean mapping quality\s*=\s*([\d,.]+)", "Mean mapping quality"),
        "GC percentage": (r"GC percentage\s*=\s*([\d,.]+)%", "% GC"),
        "general error rate": (r"general error rate\s*=\s*([\d,.]+)", "General error rate"),
        "Mean coverage": (r"mean coverageData\s*=\s*([\d,.]+)X", "Mean coverage"),
        "Std coverage": (r"std coverageData\s*=\s*([\d,.]+)X", "Std coverage")
    }

    results = []
    temp_dict = {}

    with open(file_path, "r") as fh:
        for line in fh:
            line = line.strip()

            if line.startswith("bam file ="):
                bam_path = line.split("=", 1)[1].strip()
                bam_name = os.path.basename(bam_path)
                results.append(("File name", bam_name))
                temp_dict["File name"] = bam_name

            for _, (pattern, output_name) in patterns.items():
                if output_name in dict(results):
                    continue
                match = re.search(pattern, line)
                if match:
                    value = match.group(1)
                    results.append((output_name, value))
                    temp_dict[output_name] = value

            match_dup = re.search(r"number of duplicated reads \(flagged\)\s*=\s*([\d,]+)", line)
            if match_dup:
                temp_dict["Duplicated reads (flagged)"] = match_dup.group(1)

    if "Number of reads" in temp_dict and "Duplicated reads (flagged)" in temp_dict:
        try:
            total_reads = int(temp_dict["Number of reads"].replace(",", ""))
            dup_reads = int(temp_dict["Duplicated reads (flagged)"].replace(",", ""))
            dup_pct = round(100 * dup_reads / total_reads, 2)
            results.append(("% Duplicated reads", f"{dup_pct}"))
        except Exception as e:
            print(f"[WARNING] Couldn't calculate % Duplicated reads: {e}")

    return results

def extract_pct_chimeras(picard_file):
    try:
        header = None
        with open(picard_file, "r") as pf:
            for line in pf:
                if line.startswith("CATEGORY") and ("PCT_CHIMERAS") in line:
                    header = line.strip().split('\t')
                    continue
                if line.startswith("PAIR") or line.startswith("UNPAIRED"):
                    if not header:
                        continue
                    fields = line.strip().split('\t')
                    idx = header.index("PCT_CHIMERAS")
                    value = fields[idx]
                    pct_chimeras = round(100 * float(value), 4)
                    return f"{pct_chimeras}"
    except Exception as e:
        print(f"[WARNING] Couldn't extract PCT_CHIMERAS from {picard_file}: {e}")
    return None

def extract_collect_metrics(collect_file):
    try:
        with open(collect_file, 'r') as f:
            lines = [line.strip() for line in f if not line.startswith('#') and line.strip()]
        if len(lines) >= 2:
            headers = lines[0].split('\t')
            values = lines[1].split('\t')
            data = dict(zip(headers, values))
            pf_q30 = data.get('PF_Q30_BASES', None)
            total_bases = data.get('TOTAL_BASES', None)
            return pf_q30, total_bases
    except Exception as e:
        print(f"[WARNING] Couldn't extract data from collect_bases_metrics.txt: {e}")
    return None, None

# Processing genome_results.txt
input_file =  ("genome_results.txt")
picard_file = ("picard_output.txt")

if not os.path.exists(input_file):
    print(f"[ERROR] File not found: {input_file}. Please check the input path.")
    exit(1)

parsed_results = parse_genome_results(input_file)

# Adding % Chimeras from picard_output.txt.
if os.path.exists(picard_file):
    pct_chimeras = extract_pct_chimeras(picard_file)
    if pct_chimeras:
        parsed_results.append(("% Chimeras", pct_chimeras))
    else:
        print("[INFO] % Chimeras value was not found in picard_output.txt")
else:
    print(f"[INFO] {picard_file} was not found")
  
#Adding PF_Q30_BASES and TOTAL_BASES from collect_bases_metrics.txt

collect_file = ("collect_bases_metrics.txt")

if os.path.exists(collect_file):
    pf_q30_bases, total_bases = extract_collect_metrics(collect_file)

    if pf_q30_bases:
        parsed_results.append(("PF Q30 bases", pf_q30_bases))
    else:
        print("[INFO] PF_Q30_BASES was not found in collect_bases_metrics.txt")

    if total_bases:
        parsed_results.append(("Total bases", total_bases))
    else:
        print("[INFO] TOTAL_BASES was not found in collect_bases_metrics.txt")


else:
    print(f"[INFO] {collect_file} was not found")

# Overwrite genome_results.txt file with the results. 
with open(input_file, "w") as f:
    for key, value in parsed_results:
        f.write(f"{key}\t{value}\n")


# === PART 2: Process raw_data_qualimapReport ===
raw_data_folder = ("raw_data_qualimapReport")

files_to_keep = [
    "insert_size_histogram.txt",
    "genome_fraction_coverage.txt",
    "mapped_reads_clipping_profile.txt",
    "mapping_quality_histogram.txt",
    "homopolymer_indels.txt",
    "mapped_reads_nucleotide_content.txt",
]

for filename in os.listdir(raw_data_folder):
    file_path = os.path.join(raw_data_folder, filename)
    if filename.endswith(".txt") and filename not in files_to_keep:
        try:
            os.remove(file_path)
        except Exception as e:
            print(f"[WARNING] Error al eliminar {filename}: {e}")

for filename in files_to_keep:
    original_path = os.path.join(raw_data_folder, filename)
    if not os.path.exists(original_path):
        continue

    base, ext = os.path.splitext(filename)
    if filename == "insert_size_histogram.txt":
        new_filename = "insert_size_plot_mqc.txt"
    elif filename == "mapping_quality_histogram.txt":
        new_filename = "mapping_quality_plot_mqc.txt"
    else:
        new_filename = f"{base}_mqc{ext}"

    new_file_path = os.path.join(raw_data_folder, new_filename)
    os.rename(original_path, new_file_path)

    try:
        if new_filename == "mapped_reads_nucleotide_content_mqc.txt":
            df = pd.read_csv(new_file_path, sep='\t')
            df.columns = df.columns.str.strip()
            df.rename(columns={'# Position (bp)': 'Position (bp)'}, inplace=True)

            if 'Position (bp)' in df.columns:
                for col in df.columns:
                    if col != 'Position (bp)':
                        df[col] = pd.to_numeric(df[col], errors='coerce')
                df_t = df.set_index('Position (bp)').transpose()
                df_t.to_csv(new_file_path, sep='\t', float_format='%.15f')
            else:
                print(f"[WARNING] Column 'Position (bp)' not found in {new_filename}")
        else:
            with open(new_file_path, 'r') as f:
                lines = f.readlines()

            lines = lines[1:]

            def convert_line(line):
                parts = line.strip().split()
                new_parts = []
                for part in parts:
                    try:
                        num = float(part)
                        if num.is_integer():
                            new_parts.append(str(int(num)))
                        else:
                            new_parts.append(f"{num:.15f}")
                    except ValueError:
                        new_parts.append(part)
                return '\t'.join(new_parts) + '\n'

            lines = [convert_line(line) for line in lines]

            if new_filename == "homopolymer_indels_mqc.txt":
                lines.insert(0, "Label\tNumber of Indels\n")

            with open(new_file_path, 'w') as f:
                f.writelines(lines)

    except Exception as e:
        print(f"[ERROR] Error processing {new_filename}: {e}")
print(f"[INFO] All processing complete!")
