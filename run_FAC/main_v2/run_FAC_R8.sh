#!/bin/bash
#SBATCH --job-name=FAC
#SBATCH -C rocky8
#SBATCH --array=11
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p sched_mit_psfc_r8
#SBATCH --output=Xe_v1/slurm_%a.out
#SBATCH --error=Xe_v1/slurm_%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=cjperks@mit.edu

# Init
module purge
. $HOME/bin/init_R8_py310.sh

# Your job script here
python $HOME/atomic_world/run_FAC/main_v2/run_main_v2.py $SLURM_ARRAY_TASK_ID Xe Xe_v1/

