#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cores=1
#SBATCH --time=00:05:00
#SBATCH --constraint=haswell
#SBATCH --account=m2282

module load python/3.7-anaconda-2019.07
source activate imac
srun python -u run_hmm_inference_vary_moving_window.py > run_hmm_inference_vary_moving_window.out
