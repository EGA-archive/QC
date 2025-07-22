### Execution of `BAM_pipeline_2.py` on WGS data
We executed the `BAM_pipeline_2.py` script on a whole-genome sequencing (WGS) BAM file (**175G**) from the 1000 Genomes Project dataset. If you want to download it 
click [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/high_coverage_alignment/HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam). Using a virtual machine configured with:

- **32 CPU cores**  
- **80 GB RAM**

The total execution time was as follows:

- **real: 987m29.785s**  
  Total wall-clock time (actual elapsed time from start to finish).

- **user: 984m1.076s**  
  CPU time spent in user mode (running your program's code). This can exceed real if multiple threads or CPU cores were used.

- **sys: 39m9.524s**  
  CPU time spent in system (kernel) mode, such as file I/O and other system-level operations.

👉 **In summary**: the pipeline ran for approximately **987 minutes of real time** (~16.5 hours), utilized **over 984 minutes of CPU time** through multithreading, and spent around **39 minutes in system-level tasks** like disk access and process management.
