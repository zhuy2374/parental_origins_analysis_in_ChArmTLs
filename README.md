#### ChArmTLs Parental Origin Analysis Pipeline

#### run_telogator2

This directory contains the batch processing and data aggregation scripts for running Telogator2. 

* **`batch_telogator2_jobs.sh`**
  This script iterates through a cohort sample list and dynamically generates individual Slurm submission scripts.

* **`merge_telogator2_results.sh`**
  This aggregation script collects the resulting `tlens_by_allele.tsv` files from the individual Telogator2 runs. It consolidates them into a single, unified TSV file, appending the specific sample identifiers to each row to ensure accurate tracking for subsequent statistical analysis.

#### downstream_analysis

This directory contains the R scripts used for the data cleaning, statistical analysis, and visualization. 

* **`01_qc_and_filtering.R`**
  Performs essential quality control on the raw merged data.

* **`02_post_qc_summary.R`**
  Generates descriptive statistics and visualization for the post-QC data.

* **`03_sex_age_and_ranking_conservation.R`**
  Conducts a multi-faceted statistical evaluation and visualization of ChArmTLs:
  1. Variation of ChArmTLs across different sexes.
  2. Correlation between ChArmTLs and age.
  3. Ranking conservation of ChArmTLs.

* **`04_parental_origin_comparison.R`**
  Compares the ChArmTL differences between paternal and maternal haplotypes within the offspring. This script applies statistical tests to rigorously evaluate and visualize the parent-of-origin effects on telomere lengths.

* **`05_validation_in_CEPH-1463.R`**
  Performs an independent validation analysis using the classic CEPH 1463 family dataset.
  
