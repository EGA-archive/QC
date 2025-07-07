# BAM/CRAM Quality Assessment Pipelines Documentation

## Overview
The European Genome-phenome Archive (EGA) currently stores nearly 3 million BAM and CRAM files, a number that continues to grow thanks to the tireless efforts of the scientific community. To enhance the quality reports we generate for each of these files, we've developed two custom pipelines that integrate various bioinformatics tools for comprehensive quality evaluation.

- `BAM_pipeline_2.py`: Computes a wide array of quality metrics.
- `BAM_finalize_2.py`: Post-processes the outputs and adapts them for visualization in MultiQC.

This documentation provides an in-depth overview of these tools, including the modules developed and the modifications made to integrate with MultiQC.

---

## CRAM to BAM Conversion
Working directly with CRAM files is often problematic due to inconsistent support across tools. We therefore convert CRAM files to BAM using `samtools view -b`, which requires the MD5 checksum of the reference genome unless the reference sequence is embedded in the CRAM file.

---

## Tools Used

### Qualimap
Generates comprehensive BAM statistics such as:
- Mapping quality
- Coverage
- Nucleotide composition
- Homopolymer presence

It also produces a summary file (`genome_results.txt`) with global metrics.

### RSeQC
Uses a BED annotation file to determine where reads fall (exons, introns, UTRs) and creates informative plots.

### Samtools (depth -a)
Provides per-base coverage in three columns: chromosome, position, and coverage.

### Picard (CollectAlignmentSummaryMetrics)
Generates alignment quality metrics. We extract `PCT_Chimeras` for inclusion in our summary report.

---

## Custom Functions in the Pipeline

### `parse_bam_header()`
Extracts metadata from the BAM header:
- Format version
- Sort order
- Platform
- Sample name
- Tools used

### `calculate_mapq_median()`
Calculates mapping quality (MAPQ) across the genome in 3 Mb windows. For each read, it captures:
- Chromosome
- Position
- Mapping quality

It then computes the median MAPQ per window and stores the result in a JSON file.

---

## MultiQC Modules

### Custom Tables
The default `General Stats Table` from MultiQC was insufficient. We developed a custom module that reads `genome_results.txt` and extracts relevant metrics:
- Sample name
- File name
- Total reads
- Mapped reads
- Properly paired reads
- Singleton reads
- % Duplicated reads
- Mean insert size
- Mean mapping quality
- GC content
- Error rate (mismatches)
- Mean coverage
- Coverage standard deviation
- `PCT_chimeras` from Picard

### BAM Header Table
Extracts:
- Sample name
- BAM format version
- Sort order
- Assembly ID
- Sequencing platform
- Alignment tools used

---

## Creating a New MultiQC Module

### Setup
```bash
git clone https://github.com/MultiQC/MultiQC.git
cd MultiQC
git checkout v1.25.2
git checkout -b Aurora
```

Follow the official development guide: https://docs.seqera.io/multiqc/development/modules

### Structure
```plaintext
multiqc/
└── modules/
    └── tablemaker/
        ├── __init__.py
        ├── tablemaker.py
        └── tests/
            └── test_tablemaker.py
```

### Integration
- Register the module in `pyproject.toml`
- Reinstall MultiQC locally with `pip install -e .`
- Prioritize appearance in `config_defaults.yaml`:
```yaml
module_order:
  - tablemaker
  - qualimap
```

### File Detection
Define file patterns in `search_patterns.yaml`:
```yaml
tablemaker/header:
  fn: "bam_header_info.txt"
  num_lines: 25

tablemaker/genome:
  fn: "genome_results.txt"
  num_lines: 25
```

Disable conflicting parsing by removing `genome_results.txt` from the qualimap module.

---

## Parsing Helpers
A generic function to extract key-value pairs from files:
```python
def parse_tabbed_key_value(file_contents):
    data = {}
    for line in file_contents.strip().splitlines():
        parts = line.strip().split("\t")
        if len(parts) == 2:
            key, value = parts
            data[key.strip()] = value.strip()
    return data
```

---

## Qualimap Graphs and Compatibility
Some Qualimap graphs (like mapped reads nucleotide content) are not parsed correctly by MultiQC, depending on the version. We address this in `BAM_finalize_2.py` by:
- Removing headers from `.txt` files
- Adding `Labels` and `Number of Indels` columns to `homopolymer_indels`
- Renaming files to `_mqc.txt`

---

## Genome-wide Visualization
We created two modules: `genomewide_mapq` and `genomewide_coverage`. These modules render interactive plots of mapping quality and coverage across the genome.

### Genome Coverage Calculation
We benchmarked `samtools depth -a` vs `bedtools genomecov -d`:
- Samtools is significantly faster (~10 min vs ~4 hours)
- Outputs were redirected directly to a Python script (`intervalos.py`) for 3 Mb window aggregation

### Implementation Details
- Chromosome windows are cumulative across the genome
- Chromosome lengths are extracted from `@SQ` header lines
- Chromosome names are standardized using `chrnames.json`

### Limitations
~700,000 files lack headers, making length extraction impossible. We decided not to modify the script since future submissions should include headers.

### Output
JSON with:
- Chromosome (e.g., `chr1`)
- Global window start position
- Median coverage or MAPQ

---

## Custom Graph: Median Mapping Quality
Script path: `/home/aumoreno/Desktop/bam/mediana/1000_genomes/mt.py`

Since Qualimap's internal method was unclear, we created a custom script:
- For each position, computes the median of all overlapping reads
- Aggregates in 3 Mb windows

The result replicates Qualimap’s graph and even corrects its minor inaccuracies (e.g., chromosome boundaries).

---

## How to Run MultiQC
Ensure your `bam_analysis_results` folder contains:
- `bam_header_info.txt`
- `coverage_median_3Mb_complete.json`
- `mapq_median_3Mb_complete.json`
- `genome_results.txt`
- `read_distribution.txt`
- `multiqc_config.yaml`
- `raw_data_qualimapReport/`

Inside `raw_data_qualimapReport/`, include:
- `genome_fraction_coverage_mqc.txt`
- `homopolymer_indels_mqc.txt`
- `insert_size_plot_mqc.txt`
- `mapped_reads_nucleotide_content_mqc.txt`
- `mapping_quality_plot_mqc.txt`

Run MultiQC:
```bash
multiqc . --force -c bam_analysis_results/multiqc_config.yaml
```

---

## Post-Processing Modifications
- Remove percentage toggle from `homopolymer indels`
- Hide sample count labels from certain graphs
- Remove timestamps from final report

