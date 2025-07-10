# BAM/CRAM QC Pipelines Documentation

## Overview
The European Genome-phenome Archive (EGA) currently stores nearly 3 million BAM and CRAM files, a number that continues to grow thanks to the tireless efforts of the scientific community. To enhance the quality reports we generate for each of these files, we've developed two custom pipelines that integrate various bioinformatics tools for comprehensive quality evaluation.

- `BAM_pipeline_2.py`: Computes a wide array of quality metrics.
- `BAM_finalize_2.py`: Post-processes the outputs and adapts them for visualization in MultiQC.

This documentation provides an in-depth overview of these tools, including the modules developed and the modifications made to integrate with MultiQC.

## BAM_pipeline_2.py

### Tools used

**Qualimap**: generates comprehensive BAM statistics, including coverage, mapping quality, GC content, homopolymers, etc. Also produces genome_results.txt with general metrics.

**RSeQC**: analyzes the distribution of reads across genomic regions (exons, introns, UTRs...) using a BED file.

**Samtools (depth -a)**: outputs base-level coverage in three columns: chromosome, position, and coverage.

**Picard**:

- CollectAlignmentSummaryMetrics: reports alignment quality metrics.
- CollectQualityYieldMetrics: provides quality and yield metrics for bases passing quality thresholds.

### Custom functions implemented

## `convert_cram_to_bam()`

Converts .cram files to .bam using samtools, if necessary.

## `create_bam_index()`

Generates a BAM index (.bai) if not present.

## `parse_bam_header()`
Extracts metadata from the BAM header (version, sort order, platform, sample, tools, assembly).

## `calculate_stat_from_counter()`

This is a helper function used inside `calculate_coverage_stat_streaming()` to compute a statistic (currently, the median) from a set of coverage values.

The input is a `Counter`, where each key represents a coverage value (e.g., 0, 1, 2...), and the associated value indicates how many times it appears in a genome window.  
This function transforms that dictionary into a list of values repeated according to their frequency, and calculates the median. If no values are present, it returns `0.0`.

---

## `calculate_coverage_stat_streaming()`

This function divides the genome into fixed-size intervals (e.g., 3 million base pairs).  
For each interval or window, it calculates the **median** of the coverage observed in that segment.

### 📌 How does it work?

The main idea is to **summarize BAM coverage** in fixed-size blocks, called **windows**, along the genome. In this case, each window spans 3,000,000 base pairs.

---

### 🔹 Step by step

#### 1. Read line by line with `samtools depth`

Command used:

```
samtools depth -a file.bam
```

This command returns three columns for each position in the genome:

```
chromosome   position   coverage
chr1         0          12  
chr1         1          13  
chr1         2          15  
...
```

Coverage represents the number of reads covering that position.

---

#### 2. Assign each position to a window

Each position is translated into a **global position** (adding chromosome offsets), and from that position, its window is determined.

```python
window_start = (global_pos // window_size) * window_size
```

**Example:**  
Global position `4,200,000` belongs to the window `[3,000,000 – 5,999,999]`, which starts at `3,000,000`.

---

#### 3. Accumulate the frequency of each coverage value per window

A dictionary `coverage_by_window` is built where:

- The **key** is the window start (`window_start`)
- The **value** is a `Counter`, which stores how many times each coverage level appears in that window.

**Simplified example:**

```
position   coverage
3,000,005  0
3,000,010  0
3,000,050  1
3,000,300  1
3,000,700  2
```

Then:

```python
coverage_by_window[3000000] = Counter({0: 2, 1: 2, 2: 1})
```

---

### 🔁 How is all that summarized?

Once all `Counter`s are built per window, a single representative value per window is computed: **median** coverage.

This is done with `calculate_stat_from_counter(counter)`.

```python
# Input:
Counter({0: 2, 1: 2, 2: 1})

# Converted to list:
[0, 0, 1, 1, 2]

# Median:
1
```

This value is linked to the window in the final result.

---

## 🔢 What is a global position?

The idea is to represent genome coverage as if it were a **single continuous line**. Each position must have a unique value indicating its place in this “concatenated genome”.

To achieve this, a dictionary of **accumulated offsets** is built where:

- `chr1` starts at position `0`
- `chr2` starts right after `chr1` ends (e.g., at position `249,250,621`)
- `chr3` starts after `chr2`, and so on

That way, if a window starts at global position `3,000,000`, it can be mapped to `chr1`, `chr2`, etc., depending on these offsets.  
It’s a practical way to treat the genome as a single sequence.

---

## 🔀 Chromosome name handling

A significant part of the script is dedicated to managing the various formats of chromosome names found in BAM files (e.g., `1` vs `chr1`, `MT` vs `chrM`, etc.).  

Since there’s no single standard (NCBI, UCSC, Ensembl use different conventions), a file called `chrnames.json` is used, containing a dictionary of equivalents.

If a BAM file has unknown names, the script simply uses the header names and includes them as-is in the final JSON.

---

## ⚠️ Special treatment of `chrM`

The mitochondrial chromosome (`chrM`) is very small (e.g., 16,569 bp), so it’s **not grouped into windows**.  

Instead, a **single coverage value** (median) is computed for it and placed at the end of the JSON.

---

## 🧾 Final output

The output is a `.json` file with the following format:

```json
[
  ["chr1", 0, 24.2],
  ["chr1", 3000000, 23.8],
  ["chr2", 249250621, 27.5],
  ...
  ["chrM", 900000000, 98.1]
]
```

Each line represents:

- The **chromosome** to which the window belongs  
- The **global position** in the genome (computed using offsets)  
- The **coverage value** for that window (median)

---

## BAM_finalize_2.py

### Part 1: Create modified `genome_results.txt`

Merges information from:

* `genome_results.txt` (Qualimap)
* `picard_output.txt` and `collect_bases_metrics.txt` (Picard)

To generate a new `genome_results.txt` with the following variables:

* Sample (from BAM header)
* File name 
* Number of reads (from Qualimap)
* Number of mapped reads (from Qualimap)
* Number of mapped paired reads (both in pair) (from Qualimap)
* Number of mapped paired reads (singletons) (from Qualimap)
* Median Insert size (from Qualimap)
* Mean mapping quality (from Qualimap) 
* % GC (from Qualimap) 
* General error rate (from Qualimap)
* Mean coverage (from Qualimap)
* Std coverage (from Qualimap)
* % Duplicated reads (from Qualimap) 
* PCT_chimeras (from CollectAlignmentSummaryMetrics)
* PF_Q30_BASES (from CollectQualityYieldMetrics)
* TOTAL_BASES (from CollectQualityYieldMetrics)

### Part 2: Format Qualimap files for MultiQC

**Deleted files:** removes all `.txt` files **not in**:

* `insert_size_histogram.txt`
* `genome_fraction_coverage.txt`
* `mapped_reads_clipping_profile.txt`
* `mapping_quality_histogram.txt`
* `homopolymer_indels.txt`
* `mapped_reads_nucleotide_content.txt`

**Renamed for MultiQC:**

* `insert_size_histogram.txt` → `insert_size_plot_mqc.txt`
* `mapping_quality_histogram.txt` → `mapping_quality_plot_mqc.txt`
* Others → `<basename>_mqc.txt`

**Content reformatting:**

* Removes the first line from each file.
* Converts numeric values:

  * Integers → `int`
  * Decimals → `float` with 15 decimal places
* `homopolymer_indels_mqc.txt`: adds header `Label\tNumber of Indels`
* `mapped_reads_nucleotide_content_mqc.txt`:

  * Renames `# Position (bp)` → `Position (bp)`
  * Transposes the table

---

## MultiQC Integration

### Developed/modified modules:

* `tablemaker/` (custom):

  * Parses `bam_header_info.txt` and `genome_results.txt`
  * Generates "BAM Header Info" and "Genome Results Summary" tables

* `genomewide_coverage/` (custom):

  * Visualizes genome-wide coverage across chromosomes
  * Input: `coverage_median_3Mb_complete.json`

* Modified Qualimap outputs:
  
  Graphics:
  
   * insert size distribution
   * genome fraction coverage
   * mapped reads clipping profile
   * mapping quality distribution
   * homopolymer indels
   * mapped reads nucleotide content

  * Original files parsed:

    * `insert_size_histogram.txt`
    * `genome_fraction_coverage.txt`
    * `mapped_reads_clipping_profile.txt`
    * `mapping_quality_histogram.txt`
    * `homopolymer_indels.txt`
    * `mapped_reads_nucleotide_content.txt`
  * See Qualimap docs for further explanation: [http://qualimap.conesalab.org/doc_html/analysis.html#bam-qc](http://qualimap.conesalab.org/doc_html/analysis.html#bam-qc)

* `read_distribution` (modified):

  * Input: `read_distribution.txt` (RSeQC)
  * Adds detailed descriptions based on BAM/CRAM generation strategy

### Visual configuration (`multiqc_config.yaml`):

* Configures axis labels and descriptions for the above plots
* Removes:

  * Percentage toggle button in `homopolymer indels`
  * Sample count labels in `homopolymer indels` and `mapped nucleotide content`
  * Report timestamp and input path
