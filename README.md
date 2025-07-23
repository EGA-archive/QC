# BAM/CRAM QC

The European Genome-phenome Archive (EGA) currently stores nearly 3 million BAM and CRAM files — and this number continues to grow thanks to the contributions of the scientific community. 

To improve the quality reports we generate for each of these files, we have developed a set of pipelines that automate the use of multiple bioinformatics tools for comprehensive quality assessment.

If you'd like to use these pipelines, please follow the [Start Guide](https://github.com/EGA-archive/BAM_QC/blob/main/docs/Start_Guide.md). For further details on how the scripts work, refer to the [Documentation](https://github.com/EGA-archive/BAM_QC/blob/main/docs/documentation.md).

---

## Example: BAM file from the 1000 Genomes Project

To illustrate the pipeline, we ran it on a small BAM file (Note: The BAM file contains alignments exclusively from chromosome 11) from the 1000 Genomes Project.  
You can download the input file [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam), and you can review the resulting [output folder](test/output).  
It matches the structure and content you should expect if you follow the steps in the guide.

To learn more about the 1000 Genomes Project, visit their [official website](https://www.internationalgenome.org/).

---

## Running the QC pipeline with Docker

To simplify installation and avoid dependency issues, we provide a Docker-based setup that runs the entire pipeline end-to-end.

### 1. Build the Docker image

From the root of the repository:

```bash
docker build -t bam-qc .
```

### 2. Run the pipeline

If your BAM file is located at `/absolute/path/to/muestra1.bam`, run the container like this:

```bash
docker run --rm \
  -v /absolute/path/to:/data \
  -v $(pwd)/output:/app/output \
  bam-qc /data/muestra1.bam
```

This command:

- Mounts your local folder containing the BAM file as `/data` inside the container
- Mounts the repository's `output/` directory (already created as `BAM_QC/output/`) to store the results
- Executes the full QC workflow and generates a `multiqc_report.html` in the `output/` folder

You can open the HTML report with any web browser.


---

## Benchmarking on other file types

We know the test file is relatively small, so we also evaluated the pipeline on:

- A **176 GB WGS BAM** file
- A **5.8 GB RNA-seq BAM** file

You can check the runtime performance and resource usage in the [`test/performance_logs`](test/performance_logs) folder.







