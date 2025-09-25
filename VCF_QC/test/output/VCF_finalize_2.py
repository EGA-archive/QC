#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
VCF_finalize_2.py: Generate a summary JSON report from VCF QC outputs.
"""

import os
import sys
from pathlib import Path
import pandas as pd
import json

# Global variables to collect metrics
nSamples = 0
NumRecords = 0
NumSNPs = 0
MultiAlSites = 0
MultiAlSNPSites = 0
NumInsertions = 0
NumDeletions = 0
Ts = 0
Tv = 0
TsTvRatio = 0
IndelDist = {}
BaseChangeDist = {"A": {}, "T": {}, "C": {}, "G": {}}
QualityDist = {}
AverageQuality = 0
Version = "Na"
DepthDist = {}
PSC = {}
PSI = {}
stats_table = {}
FrequencyDist = {}
AverageDepth = 0


def parse_stats(lines):
    """Parse bcftools stats output to extract variant statistics."""
    global nSamples, NumRecords, NumSNPs, MultiAlSites, MultiAlSNPSites
    global NumInsertions, NumDeletions, Ts, Tv, TsTvRatio, IndelDist
    global BaseChangeDist, QualityDist, AverageQuality
    global FrequencyDist, DepthDist, AverageDepth
    global PSC, PSI, stats_table

    quality_data, af_data, depth_data, psc_data, psi_data = [], [], [], [], []
    BaseChangeDist = {"A": {}, "T": {}, "C": {}, "G": {}}
    IndelDist = {}
    QualityDist = {}
    FrequencyDist = {}
    DepthDist = {}

    for line in lines:
        fields = line.strip().split("\t")
        if not fields or len(fields) < 2:
            continue

        if fields[0] == "SN":
            if fields[2] == "number of samples:":
                nSamples = int(fields[3])
            elif fields[2] == "number of records:":
                NumRecords = int(fields[3])
            elif fields[2] == "number of SNPs:":
                NumSNPs = int(fields[3])
            elif fields[2] == "number of multiallelic sites:":
                MultiAlSites = int(fields[3])
            elif fields[2] == "number of multiallelic SNP sites:":
                MultiAlSNPSites = int(fields[3])

        elif fields[0] == "IDD":
            val = int(fields[2])
            count = int(fields[3])
            if val < 0:
                NumDeletions += count
            elif val > 0:
                NumInsertions += count
            IndelDist[val] = count

        elif fields[0] == "TSTV":
            Ts = int(fields[2])
            Tv = int(fields[3])
            TsTvRatio = float(fields[4])

        elif fields[0] == "ST":
            ref, alt = fields[2].split(">")
            BaseChangeDist[ref][alt] = int(fields[3])

        elif fields[0] == "QUAL" and fields[2] != ".":
            q = float(fields[2])
            count = int(fields[3])
            quality_data.append([q, count])
            QualityDist[q] = count

        elif fields[0] == "AF" and nSamples >= 50:
            af = float(fields[2])
            count = int(fields[3])
            af_data.append([af, count])
            FrequencyDist[af] = count

        elif fields[0] == "DP":
            dp_str = fields[2].strip()
            dp = 500 if dp_str == ">500" else float(dp_str)
            count = int(fields[3])
            if count > 0:
                depth_data.append([dp, count])
                DepthDist[dp] = count

        elif fields[0] == "PSC":
            sample = fields[2]
            n_ref_hom = int(fields[3])
            n_non_ref_hom = int(fields[4])
            n_hets = int(fields[5])
            n_transitions = int(fields[6])
            n_transversions = int(fields[7])
            avg_depth = float(fields[9])
            n_missing = int(fields[1]) if nSamples >= 50 else "null"
            ratio_het_hom = n_hets / (n_ref_hom + n_non_ref_hom) if (n_ref_hom + n_non_ref_hom) > 0 else None
            ratio_ts_tv = n_transitions / n_transversions if n_transversions > 0 else None
            psc_data.append([sample, avg_depth, n_missing, ratio_het_hom, ratio_ts_tv])

        elif fields[0] == "PSI":
            sample = fields[2]
            psi_data.append([
                sample,
                int(fields[7]),
                int(fields[8]),
                int(fields[9]),
                int(fields[10]),
            ])

    PSC = pd.DataFrame(psc_data, columns=["sample", "average_depth", "n_missing", "ratio_het_hom", "ratio_ts_tv"])
    PSI = pd.DataFrame(psi_data, columns=["sample", "nInsHets", "nDelHets", "nInsAltHoms", "nDelAltHoms"])
    stats_table = pd.merge(PSC, PSI, on="sample") if not PSC.empty and not PSI.empty else pd.DataFrame()

    QualityDist_df = pd.DataFrame(quality_data, columns=["Quality", "Number of SNPs"])
    DepthDist_df = pd.DataFrame(depth_data, columns=["Depth", "Number of SNPs"])
    AverageQuality = QualityDist_df["Quality"].mean() if not QualityDist_df.empty else 0
    AverageDepth = DepthDist_df["Depth"].mean() if not DepthDist_df.empty else 0


def make_json():
    """Create the JSON summary from collected metrics."""
    indels = [{"Indel size": k, "Number of SNPs": v} for k, v in sorted(IndelDist.items())]
    quality = [{"Quality": k, "Number of SNPs": v} for k, v in sorted(QualityDist.items())]
    depth = [{"Depth": k, "Number of SNPs": v} for k, v in sorted(DepthDist.items())]
    freq_agg = {}
    for k, v in FrequencyDist.items():
        key = int(k * 100)
        freq_agg[key] = freq_agg.get(key, 0) + v
    freq = [{"Frequency (%)": k, "Number of SNPs": v} for k, v in sorted(freq_agg.items())]
    stats = [
        {k: round(v, 2) if isinstance(v, (int, float)) else v for k, v in row.items()}
        for row in stats_table.to_dict(orient="records")
    ] if not stats_table.empty else []

    return json.dumps({
        "Report Version": 2,
        "VCF Version": Version,
        "NbSamples": nSamples,
        "NbRecords": NumRecords,
        "AvgQuality": AverageQuality,
        "MultiAlSites": MultiAlSites,
        "AvgDepth": AverageDepth,
        "MultiAlSNPSites": MultiAlSNPSites,
        "NumSNPs": NumSNPs,
        "Variants": [
            {"SNPs": NumSNPs},
            {"Ins": NumInsertions},
            {"Del": NumDeletions},
        ],
        "TsTvRatio": TsTvRatio,
        "QualityDist": quality,
        "Depth": depth,
        "Stats_table": stats,
        "FrequencyDist": freq,
        "IndelDistribution": indels
    }, indent=4)


def main():
    wd = Path.cwd()
    output_path = wd.parent / "output" / "final_report.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    header_path = wd / "header.txt"
    version_path = wd / "vcfVersion.txt"
    stats_sample = wd / "stats_sample.txt"
    stats_default = wd / "stats.txt"

    if not header_path.exists():
        print("[ERROR] header.txt not found.")
        return

    stats_path = stats_sample if stats_sample.exists() and stats_sample.stat().st_size > 0 else stats_default
    if not stats_path.exists() or stats_path.stat().st_size == 0:
        print("[ERROR] No stats_sample.txt or stats.txt found.")
        return

    with open(stats_path) as f:
        parse_stats(f)

    global Version
    if version_path.exists():
        Version = version_path.read_text().strip()

    report = make_json()
    with open(output_path, 'w') as f:
        f.write(report)

    print(f"[INFO] Final report written to {output_path}")


if __name__ == "__main__":
    assert sys.version_info >= (3, 6), "This script requires Python >= 3.6"
    main()
