# Quick Start Guide: BAM/CRAM Quality Assessment
 
## 1. Requirements

Before starting, make sure you have the following:
# Files

* A BAM or CRAM file (of course!)
* A FASTA reference file that matches the reference assembly of the BAM/CRAM
* A BED file with annotations (same assembly as the BAM/CRAM)
* `chrnames.json` (download from the `Required_Files` folder in this repository).
  This must be placed in the same directory where you are running the BAM_pipeline_2.py script.
* `multiqc_config.yaml` (also from `Required_Files`). Needed for generating the MultiQC report.


# Programming language 
* python3
  
# Tools
* MultiQC-EGA (custom fork of MultiQC adapted for this pipeline).
* RSeQC
* Qualimap
* samtools
* Picard
  
# Libraries
* subprocess
* argparse
* os
* pysam
* json
* math
* statistics
* collections
* re
* pandas
* warnings

---

## 2. Installing MultiQC-EGA

From your terminal, clone the customized MultiQC repository:

```bash
git clone https://github.com/EGA-archive/BAM_QC.git
cd MultiQC-EGA
pip install -e .
```
You're now ready to use MultiQC-EGA.

---

## 3. Running the Pipeline

### Step 1: Run `BAM_pipeline_2.py`

It is located in the scripts folder in this repository. Remember to place chrnames.json file where you are running this script. 
```bash
./BAM_pipeline_2.py --bam <your_input.bam_or_cram> --bed <annotations.bed> --fasta <reference.fa>
```

This script will:

* Convert CRAM to BAM (if needed)
* Extract header metadata
* Run QC tools: Qualimap, Picard, samtools, read\_distribution
* Compute MAPQ and coverage stats
* Generate results in a folder named `bam_analysis_results/` (same directory as your input file)

---

### Step 2: Post-process with `BAM_finalize_2.py`

Run it in the same folder where `bam_analysis_results/` was generated:

```bash
./BAM_finalize_2.py
```

> This script modifies `genome_results.txt` and files inside `raw_data_qualimapReport/`.
>
> If you'd like to preserve the original Qualimap files, create a backup folder (e.g. `Qualimap/`) and move them there before running this step.

---

### Step 3: Generate the MultiQC Report

Make sure your working directory contains:

* `bam_analysis_results/` (with all generated files inside)
* `multiqc_config.yaml` (in the same directory, not inside `bam_analysis_results/`)

Then run:

```bash
multiqc . -c multiqc_config.yaml -e picard -e qualimap --force
```

This will generate a customized, interactive MultiQC report summarizing all your QC results.

---

## Final Notes

* The `chrnames.json` file ensures chromosome naming is standardized across tools.
* If anything fails, check that your reference genome and annotation files match the same genome assembly as the BAM/CRAM.
