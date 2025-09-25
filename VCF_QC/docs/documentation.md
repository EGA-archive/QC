This repository provides a two-step Python pipeline to extract and summarize quality control (QC) metrics from VCF or VCF.GZ files. The pipeline is designed for reproducible and interpretable processing of genomic variant data.

---

## Pipeline Components

### 1. `run/VCF_pipeline_2.py`

**Main Steps**:

1. **Input validation**: Accepts `.vcf` or `.vcf.gz` as input.
2. **Indexing**: Creates a `.csi` index using `bcftools`.
3. **Random SNP selection**: Extracts 1000 random `rsID` SNPs.
4. **Sample listing**: Generates a list of samples using `bcftools query -l`.
5. **Chromosome list**: Extracts unique chromosome names.
6. **VCF version detection**: Parses the first line of the VCF header.
7. **Header extraction**: Saves all header lines excluding `#CHROM`.
8. **Sample-level statistics**: Uses `bcftools stats -S` to extract per-sample metrics.
9. **General statistics**: Runs `bcftools stats` for full VCF-wide metrics.

**Outputs** (saved in `output/`):

- `<sample>.csi`: VCF index
- `random_snps.txt`: 1000 random rsID SNPs
- `samples.txt`: List of samples in the VCF
- `chromosomes.txt`: Unique chromosomes
- `version.txt`: Detected VCF format version
- `header.txt`: Full VCF header excluding #CHROM
- `stats_sample.txt`: Per-sample statistics (if available)
- `stats.txt`: Overall file-level variant statistics

---

### 2. `output/VCF_finalize_2.py`

**Purpose**: Parse and combine intermediate files into a single structured JSON report, `final_report.json`.

If the `stats_sample.txt` file is empty (typically due to a missing `sample_list.txt` file), the script will instead parse `stats.txt`, and the `PerSampleStats` section (explained below) will be left empty.


## `final_report.json` – Field Descriptions

The `final_report.json` file contains structured quality metrics extracted from the input VCF file using `bcftools`. It is organized into several sections:

---

### PerSampleStats

A list of entries (one per sample) containing summary statistics:

| Field             | Description |
|-------------------|-------------|
| `sample`          | Sample identifier. |
| `average_depth`   | Average sequencing depth across variants (from FORMAT/DP). |
| `n_missing`       | Number of missing genotypes (`./.`) in the sample. |
| `ratio_het_hom`   | Ratio of heterozygous to homozygous alternate genotypes. |
| `ratio_ts_tv`     | Transition/transversion ratio for the sample. |
| `nInsHets`        | Number of heterozygous insertions. |
| `nDelHets`        | Number of heterozygous deletions. |
| `nInsAltHoms`     | Number of homozygous alternate insertions. |
| `nDelAltHoms`     | Number of homozygous alternate deletions. |

---

### FrequencyDist

Distribution of alternate allele frequencies across all variants:

| Field             | Description |
|-------------------|-------------|
| `Frequency (%)`   | Rounded alternate allele frequency (in %). |
| `Number of SNPs`  | Number of SNPs at that frequency level. |

This can be used to generate allele frequency histograms.

---

### QualityDist

Distribution of variant quality scores (`QUAL` field):

| Field   | Description |
|---------|-------------|
| `QUAL`  | Phred-scaled variant quality score. |
| `count` | Number of variants with that score. |

Useful for deciding quality thresholds for filtering.

---

### DepthDist

Distribution of total sequencing depth (`DP`) across variants:

| Field | Description |
|-------|-------------|
| `DP`  | Total depth of coverage for the variant. |
| `count` | Number of variants at that depth. |

Helps identify under- or over-covered regions.

---

### IndelTypeDistribution

Counts of insertion and deletion variants:

| Field    | Description |
|----------|-------------|
| `Type`   | Either `Insertion` or `Deletion`. |
| `Count`  | Total number of indels of that type. |

Can be used to check for balance or bias in indel calling.

---

### GeneralStats

Overall summary metrics for the VCF file:

| Field                      | Description |
|----------------------------|-------------|
| `Number of records`        | Total number of variant entries in the file. |
| `Number of samples`        | Total number of samples. |
| `Number of SNPs`           | Total number of single-nucleotide polymorphisms. |
| `Number of indels`         | Total number of insertions and deletions. |
| `Number of transitions`    | Count of transition mutations (A↔G, C↔T). |
| `Number of transversions`  | Count of transversion mutations. |
| `ts/tv ratio`              | Transition/transversion ratio. |
| `Percent missing genotypes` | Percentage of genotypes that are missing. |

These values give a global snapshot of the VCF file's quality and content.

---


## Requirements

- Python ≥ 3.6
- Tools: `bcftools`, `awk`, `grep`, `shuf`, `zcat`, `cat`
- Python packages: `pandas`, `json`, `argparse`
