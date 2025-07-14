# Quick Start Guide: BAM/CRAM Quality Assessment
 
## 1. Requirements

Before starting, make sure you have the following:
# Files

* A BAM or CRAM file (of course!)
* A FASTA reference file that matches the reference assembly of the BAM/CRAM
* A BED file with annotations (same assembly as the BAM/CRAM)
* `chrnames.json` (you already have it in place).
* `multiqc_config.yaml` (also in placed). Needed for generating the MultiQC report.


# Programming language 
* python3


---

## 2. Installing required tools and libraries

### Step 1: Download this repository to your computer.

https://docs.github.com/es/get-started/start-your-journey/downloading-files-from-github

### Step 2: Enter run folder

Now you can start working on your terminal, entering the run folder. 

```bash
cd run/
```

### Step 3: Install libraries, RSeQC and MULTIQC-EGA 

```bash
pip install -r requirements.txt
```
You're now half way to start running our pipeline!

---
### Step 4: Install samtools 

This is not a python tool, so we could not include it in the requirements.txt file, and we could not include it in this repository wither, so you need to manage it by yourself, but it should be easy, you are almost there.
Follow the tool documentation: https://github.com/samtools/samtools?tab=readme-ov-file

You will end up having: 

# Tools
* MultiQC-EGA (custom fork of MultiQC adapted for this pipeline, included in the requirements.txt).
* RSeQC (included in requirements.txt)
* Qualimap (included in this repository)
* samtools (downloaded by the user)
* Picard (included inside Qualimap folder in this repository)
  
# Libraries (included in requirements.txt)
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
  
## 3. Running the Pipeline

### Step 1: Run `BAM_pipeline_2.py`

It is located in the run folder in this repository. 

```bash
./BAM_pipeline_2.py --bam <your_input.bam_or_cram> --bed <annotations.bed> --fasta <reference.fa>
```

This script will:

* Convert CRAM to BAM (if needed)
* Create bam index file
* Extract header metadata
* Run QC tools: Qualimap, Picard, samtools, read_distribution
* Compute coverage stats
* Generate results in the folder named `output/` 

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
