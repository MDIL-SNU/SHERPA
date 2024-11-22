#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --job-name="example"


mpirun -np $SLURM_NTASKS /path/to/sherpa_ase
