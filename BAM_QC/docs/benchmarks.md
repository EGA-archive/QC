The script BAM_pipeline_2.py can take a lot of resources, if you want to learn more about how many resources took for us in different BAM files, keep reading. 

We are executing this pipeline in an Ubuntu 20.4 VM with a 32 cores CPU and 80gb of RAM.

### RNAseq file: 

We tested this [file](https://www.ega-archive.org/datasets/EGAD00001002256), 

Total time (actual time elapsed).
→ 152m48.964s: the full duration of the process.

user: CPU time spent in user mode (your program's own code).
→ 200m44.345s: can be greater than real if multiple CPU cores were used.

sys: CPU time spent in system (kernel) mode (e.g., file I/O).
→ 4m45.965s: indicates time spent in system-level operations.

In short: the process took 152 minutes of real time, used over 200 minutes of CPU time (due to multithreaded), and spent little time in system calls.

### WGS file: 

