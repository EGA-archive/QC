# BAM/CRAM QC

The European Genome-phenome Archive (EGA) currently stores nearly 3 million BAM and CRAM files — and this number continues to grow thanks to the contributions of the scientific community. 

To improve the quality reports we generate for each of these files, we have developed a set of pipelines that automate the use of multiple bioinformatics tools for comprehensive quality assessment.

If you'd like to use these pipelines, please follow the [Start Guide](https://github.com/EGA-archive/BAM_QC/blob/main/docs/Start_Guide.md). For further details on how the scripts work, refer to the [Documentation](https://github.com/EGA-archive/BAM_QC/blob/main/docs/documentation.md).

---

### Example: BAM file from the 1000 Genomes Project

To illustrate the pipeline, we ran it on a small BAM file from the 1000 Genomes Project.  
You can download the input file [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam), and you can review the resulting [output folder](test/output).  
It matches the structure and content you should expect if you follow the steps in the guide.

To learn more about the 1000 Genomes Project, visit their [official website](https://www.internationalgenome.org/).

---

### Benchmarking on other file types

We know the test file is relatively small, so we also evaluated the pipeline on:

- A **176 GB WGS BAM** file
- A **5.8 GB RNA-seq BAM** file

You can check the runtime performance and resource usage in the [`test/performance_logs`](test/performance_logs) folder.






