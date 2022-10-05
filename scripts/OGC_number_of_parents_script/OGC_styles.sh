#!/bin/bash
#SBATCH -p main
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 0-02:00:00
#SBATCH -J OGC_styles_simple_investigation
#SBATCH -o OGC_styles_simple_investigation.out
#SBATCH --array=1-250
#SBATCH --output=slurm-%A_%a.out
#SBATCH -e OGC_styles_simple_investigation.err
#SBATCH --mail-type=ALL

### Initialization
# Get Array ID
i=${SLURM_ARRAY_TASK_ID}
NAME_PARAMETER_FILE="OGC_styles_simple_investigation1"

srun hostname
sstat -j $SLURM_JOB_ID

# Loading modules
module unload gcc
module load gcc/11.2.0
module load R/4.1.2
#module load intel/oneapi/mkl

Rscript /home/WUR/nieho006/Optimum_contribution/OGC_styles_simple_investigation.R $i $NAME_PARAMETER_FILE


echo $SLURM_ARRAY_TASK_ID
echo '******************** FINISHED ***********************'
echo