# BAM/CRAM QC Pipelines Documentation

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


