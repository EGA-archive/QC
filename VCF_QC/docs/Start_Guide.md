
# VCF QC Pipeline – Getting Started Guide

This project provides a Python-based pipeline to extract quality control (QC) metrics from VCF and BCF (both uncompressed and BGZF-compressed) files and consolidate them into a final JSON report. It consists of two scripts:

- `run/VCF_pipeline_2.py`: performs preprocessing of the input VCF and generates multiple intermediate files in the `output/` folder.
- `output/VCF_finalize_2.py`: parses the outputs and produces a structured summary file named `final_report.json`.

---


## Requirements

- Python ≥ 3.6
- Dependencies:
  - `bcftools` (available in the system PATH)
  - `awk`, `grep`, `shuf`, `zcat` or `cat`
  - Python packages: `pandas`, `json`, `argparse`, etc.

Install Python dependencies with:

```bash
pip install pandas
```

Install `bcftools` via conda or apt:

```bash
conda install -c bioconda bcftools
# or
sudo apt install bcftools
```

---

## Step 1 – Run `VCF_pipeline_2.py`

This script generates QC metrics from a VCF (or VCF.GZ) file.

```bash
cd run/
./VCF_pipeline_2.py --vcf path/to/input.vcf.gz
```

 **Output files in `output/`**:

- `*.csi`: VCF index
- `random_snps.txt`: 1000 randomly selected rsID SNPs
- `samples.txt`: list of samples
- `chromosomes.txt`: list of unique chromosomes
- `version.txt`: VCF version (format)
- `header.txt`: VCF header (excluding #CHROM)
- `stats_sample.txt`: per-sample statistics (if applicable)
- `stats.txt`: global VCF statistics

---

## Step 2 – Run `VCF_finalize_2.py`

This script consolidates the previous results into a comprehensive JSON report.

```bash
cd output/
./VCF_finalize_2.py
```

This will generate:

```
output/final_report.json
```

The resulting `final_report.json` is designed for easy visualization and downstream analysis. It includes structured data on:

- **Per-sample summary statistics**, such as:
  - Average sequencing depth
  - Number of missing genotypes
  - Transition/transversion ratio (ts/tv)
  - Heterozygous/homozygous ratio
  - Indel counts split by zygosity (heterozygous/homozygous insertions and deletions)

- **Depth distribution** 
- **Quality score distribution** 
- **Allele frequency distribution**
- **Indel type distribution**

This JSON file is suitable for generating summary tables and visual dashboards. At EGA, we use the same structure to create visual QC reports at scale. If you’d like to generate similar plots or need help interpreting the output, feel free to contact us — we’d be happy to assist.
