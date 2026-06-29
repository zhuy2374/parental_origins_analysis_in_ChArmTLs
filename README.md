#### ChArmTLs Parental Origin Analysis Pipeline

#### run_telogator2

This directory contains the batch processing and data aggregation scripts for running Telogator2. 

* **`batch_telogator2_jobs.sh`**
  This script iterates through a cohort sample list and dynamically generates individual Slurm submission scripts.

* **`merge_telogator2_results.sh`**
  This aggregation script collects the resulting `tlens_by_allele.tsv` files from the individual Telogator2 runs. It consolidates them into a single, unified TSV file, appending the specific sample identifiers to each row to ensure accurate tracking for subsequent statistical analysis.
