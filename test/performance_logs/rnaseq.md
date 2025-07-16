### Execution of `BAM_pipeline_2.py` on RNA-seq data

We executed the `BAM_pipeline_2.py` script on an RNA-seq file (5.8G) from the dataset [EGAD00001002256](https://www.ega-archive.org/datasets/EGAD00001002256), using a virtual machine configured with:

- **32 CPU cores**
- **80 GB RAM**

The total execution time was as follows:

- **real**: `152m48.964s`  
  Total wall-clock time (actual elapsed time from start to finish).

- **user**: `200m44.345s`  
  CPU time spent in user mode (running your program's code). This can exceed `real` if multiple threads or CPU cores were used.

- **sys**: `4m45.965s`  
  CPU time spent in system (kernel) mode, such as file I/O and other system-level operations.

👉 **In summary**: the pipeline ran for approximately **152 minutes of real time**, utilized over **200 minutes of CPU time** thanks to multithreading, and spent relatively little time in system-level operations.
