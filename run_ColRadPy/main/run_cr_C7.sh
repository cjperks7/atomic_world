#!/bin/bash
#SBATCH --job-name=ColRadPy
#SBATCH -C centos7
#SBATCH --array=1-3
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p sched_mit_psfc
#SBATCH --output=Kr_v1/slurm_%a.out
#SBATCH --error=Kr_v1/slurm_%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=cjperks@mit.edu

# Init
module purge
. $HOME/bin/init_C7_py39.sh

# Your job script here
python $HOME/atomic_world/run_ColRadPy/main/sim_FAC_PECs_HPC.py $SLURM_ARRAY_TASK_ID Kr Kr_v1/

