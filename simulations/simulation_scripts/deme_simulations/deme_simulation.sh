#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cores=1
#SBATCH --time=48:00:00
#SBATCH --qos=premium
#SBATCH --constraint=haswell
#SBATCH --account=m2282

module load python/3.7-anaconda-2019.07
source activate imac
srun python -u deme_simulation.py > deme_simulation.out
