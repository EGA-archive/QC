# Quick Start Guide: BAM/CRAM Quality Assessment
 
## 1. Requirements

Before starting, make sure you have the following:
### Files

* A BAM or CRAM file (of course!)
* A FASTA reference file that matches the reference assembly of the BAM/CRAM
* A BED file with annotations (same assembly as the BAM/CRAM)
* `chrnames.json` (you already have it in place).
* `multiqc_config.yaml` (also in placed). Needed for generating the MultiQC report.


### Programming language 
* python3

---

## 2. Installing required tools and libraries

### Step 1: Download this repository to your computer.

If you do not know how to, you can click [here](https://docs.github.com/es/get-started/start-your-journey/downloading-files-from-github) and check GitHub documentation. 


### Step 2: Install libraries, RSeQC and MULTIQC-EGA 

```bash
cd BAM_QC
pip install -r requirements.txt
```
You're now half way to start running our pipeline!

### Step 3: Enter run folder

Now you can start working on your terminal, entering the run folder. 

```bash
cd run/
```
### Step 4: Install samtools 

This is not a python tool, so we could not include it in the requirements.txt file, and we could not include it in this repository wither, so you need to manage it by yourself, but it should be easy, you are almost there.
Follow the [tool documentation](https://github.com/samtools/samtools?tab=readme-ov-file).

### Step 5: Install picard

Click on this [link](https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar) and the download will start. If you already have it, and you do not want to download it again, please place the picard.jar file in this folder. If you do not want to, you can modify BAM_pipeline_2.py, with the relative path of your picard.jar file. 

Eventually, You will end up having: 

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
---
## 3. Running the Pipeline

### Step 1: Run `BAM_pipeline_2.py`

It is located in `run/` in this repository. 
You may need to give permission to the file to be executed. If that's the case, run: 

```bash
chmod +x BAM_pipeline_2.py
```
Then give permission to Qualimap by running: 

```bash
chmod +x qualimap_v2.3/qualimap
```
And then you can run the pipeline: 
```bash
./BAM_pipeline_2.py --bam <your_input.bam_or_cram> --bed <annotations.bed> --fasta <reference.fa>
```
Attention: This pipeline will set your 
This script will:

* Convert CRAM to BAM (if needed)
* Create bam index file
* Extract header metadata
* Run QC tools: Qualimap, Picard, samtools, read_distribution
* Compute coverage stats
* Generate results in the folder named `output/` 

---

### Step 2: Post-process with `BAM_finalize_2.py`

Access `output/`:

```bash
cd ../output
```
Execute the `BAM_finalize_2.py` pipeline, you may need to give permission to the file as before. This step will permanently modify the Qualimap output: 

```bash
./BAM_finalize_2.py
```

> This script modifies `genome_results.txt`, `picard_output.txt`, `collect_bases_metrics.txt` and files inside `raw_data_qualimapReport/`.
>
> If you'd like to preserve the original Qualimap files, create a backup folder (e.g. `Qualimap/`) and move them there before running this step.

---

### Step 3: Generate the MultiQC Report

In the same folder as the previous step, `output/`, run:

```bash
multiqc . -c multiqc_config.yaml -e picard -e qualimap --force
```

Congratulations, you finished! This will generate a customized, interactive MultiQC report summarizing all your QC results `multiqc_report.html`.

---

## Final Notes

* The `chrnames.json` file ensures chromosome naming is standardized across tools.
* If anything fails, check that your reference genome and annotation files match the same genome assembly as the BAM/CRAM.
