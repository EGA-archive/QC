# BAM/CRAM QC

The EGA currently stores nearly 3 million BAM and CRAM files — and this number keeps growing thanks to the collective efforts of the scientific community. To improve the quality reports we generate for each of these files, we’ve developed a set of pipelines that automates the use of several bioinformatics tools for comprehensive quality assessment.

BAM_pipeline_2.py:

Converts CRAM to BAM if needed.
Extracts metadata from the BAM header.
Computes quality metrics using tools like Qualimap, Picard, samtools, and RSeQC.
Calculates per-window statistics (MAPQ, coverage) using custom Python scripts.


BAM_finalize_2.py: 

Formats the results from BAM_pipeline_2.py to be compatible with MultiQC, using a customized version of MultiQC. You’ll also find this modified version included in this repository.


