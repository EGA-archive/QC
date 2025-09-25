# Quality Control Pipelines at EGA

The **European Genome-phenome Archive (EGA)** provides automated quality control (QC) pipelines for common genomic file formats to ensure files are complete, well-structured, and ready for downstream use.

- The EGA currently stores **a large number of phenotypic and genomic files** — and it continues to grow thanks to the contributions of the scientific community.
- This repository includes **Python-based pipelines for BAM/CRAM and VCF files**, designed to extract key quality metrics and generate useful summaries.

---

## 📂 Available QC Pipelines

- **BAM/CRAM QC** → Quality assessment for aligned read files (BAM/CRAM)
- **VCF QC** → Quality assessment for variant call files (VCF/VCF.GZ)

Each pipeline comes with example datasets, performance notes, and documentation.

---

## 🧬 BAM/CRAM QC

To improve the quality reports we generate for each of these files, we have developed a set of pipelines that automate the use of multiple bioinformatics tools for comprehensive quality assessment.

If you'd like to use these pipelines, please follow the **[Start Guide](https://github.com/EGA-archive/BAM_QC/blob/main/docs/Start_Guide.md)**. For further details on how the scripts work, refer to the **[Documentation](https://github.com/EGA-archive/BAM_QC/blob/main/docs/documentation.md)**.

---

### Example: BAM file from the 1000 Genomes Project

To illustrate the pipeline, we ran it on a small BAM file (**Note:** the BAM file contains alignments exclusively from chromosome 11) from the 1000 Genomes Project.

- Download the input file **[here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam)**  
- Review the resulting **[output folder](test/output)**

It matches the structure and content you should expect if you follow the steps in the guide.

Learn more about the 1000 Genomes Project on their **[official website](https://www.internationalgenome.org/)**.

---

### Running the BAM QC pipeline with Docker

#### 1) Build the Docker image

From the root of the repository:

```bash
docker build -t bam-qc .
```

> **Important:** Make sure the following files have execution permissions *before building* the image:
>
> - `run/qualimap_v2.3/qualimap`
> - `run/BAM_pipeline_2.py`
> - `output/BAM_finalize_2.py`
> - `run/wrapper.py`
>
> If not, set them manually using:
>
> ```bash
> chmod +x run/qualimap_v2.3/qualimap run/BAM_pipeline_2.py output/BAM_finalize_2.py run/wrapper.py
> ```

#### 2) Run the pipeline

If your BAM, BED and FASTA files are located in `/absolute/path/to/`, run the container like this:

```bash
docker run --rm   -v /absolute/path/to:/data   -v $(pwd)/output:/app/output   bam-qc   --bam /data/muestra1.bam   --bed /data/regions.bed   --fasta /data/reference.fasta
```

This command:

- Mounts your local folder containing the input files as `/data` inside the container
- Mounts the repository's `output/` directory (already created as `BAM_QC/output/`) to store the results
- Passes the required arguments (`--bam`, `--bed`, `--fasta`) to the pipeline
- Executes the full QC workflow and generates a `multiqc_report.html` in the `output/` folder

#### 3) Run MultiQC

To create the report please run inside the `output/` folder:

```bash
multiqc . -e picard -e qualimap -c multiqc_config.yaml
```

You're done! Check the results in the `multiqc_report.html` file.

---

### Benchmarking on other file types

We know the test file is relatively small, so we also evaluated the pipeline on:

- A **176 GB WGS BAM** file
- A **5.8 GB RNA-seq BAM** file

You can check the runtime performance and resource usage in the [`test/performance_logs`](test/performance_logs) folder.

---

## 🧩 VCF QC

The European Genome-phenome Archive (EGA) stores thousands of VCF files, which encode genetic variation data across samples and studies.

To ensure these files meet expected quality and formatting standards, we have developed a **Python-based pipeline** that performs initial QC checks and extracts summary metrics from **VCF and VCF.GZ** files.

If you'd like to use this pipeline, please follow the **[Start Guide](https://github.com/EGA-archive/QC/blob/main/VCF_QC/docs/Start_Guide.md)**. For detailed explanations of each step, refer to the **[Documentation](https://github.com/EGA-archive/QC/blob/main/VCF_QC/docs/documentation.md)**.

---

### Example: VCF QC on 1000 Genomes file

We also ran the VCF QC pipeline on a file from the 1000 Genomes Project:

- **[Download input VCF](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz)**
- **[View output folder](https://github.com/EGA-archive/QC/tree/main/VCF_QC/test/output)**

Runtime and memory usage details are included in the corresponding performance log **[here](VCF_QC/test/performance_log.md)**.

🔗 **[VCF file format specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf)**

---

Together, these pipelines provide a consistent framework for checking the integrity and usability of genomic data hosted at EGA.
