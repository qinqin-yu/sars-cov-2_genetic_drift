#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cores=16
#SBATCH --time=48:00:00
#SBATCH --qos=regular
#SBATCH --constraint=haswell
#SBATCH --account=m2282

module load python/3.7-anaconda-2019.07
source activate imac
srun python -u run_hmm_inference_1.py > run_hmm_inference_1.out
