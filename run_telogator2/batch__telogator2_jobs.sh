#!/bin/bash

# ==============================================================================
# Description: Batch generate and submit Telogator2 jobs.
# ==============================================================================

# ------------------ Paths and Environment (Users must modify) ------------------
# 1. Set the root working directory
WORK_DIR="/path/to/your_workdir"

# 2. Set the absolute path to the Conda initialization script
CONDA_SH_PATH="/path/to/miniconda3/etc/profile.d/conda.sh"

# 3. Set the absolute path to the Telogator2 python script
TELOGATOR_SCRIPT="/path/to/telogator2.py"

# 4. Set Slurm job scheduling parameters for the HPC cluster
SLURM_PARTITION="your_partition"  # Modify to match your target partition
# -------------------------------------------------------------------------------

# Path to the sample list file
SAMPLE_LIST="${WORK_DIR}/sample_name_list_file"

# Mapping file containing the absolute paths to the corresponding HiFi/CCS BAM files
BAM_PATH_LIST="${WORK_DIR}/sample_ccsbam_path_list_file"

# Ensure the script and log directory exist before generating scripts
mkdir -p ${WORK_DIR}/scripts
mkdir -p ${WORK_DIR}/logs


# Iterate through each individual sample in the list
for individual in `cat ${SAMPLE_LIST}` 
do
# Get BAM file path for the current individual
individual_bam_path=$(grep ${individual} ${BAM_PATH_LIST} | paste -s -d ' ')

# Generate Slurm job script
echo "#!/bin/bash
#SBATCH -J ${individual}_telogator2
#SBATCH -o ${WORK_DIR}/logs/${individual}_telogator2.log
#SBATCH -p <your_partition>
#SBATCH -N 1
#SBATCH --cpus-per-task=8

# Initialize and activate the Conda environment
source ${CONDA_SH_PATH}
conda activate telogator2

cd ${WORK_DIR}

# Run Telogator2 on HiFi data
python ${TELOGATOR_SCRIPT} -i ${individual_bam_path} -r hifi -n 2 -p 8 -o ${individual}/

" > ${WORK_DIR}/scripts/${individual}_telogator2

done


