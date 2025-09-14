#!/bin/bash
#SBATCH --job-name=Dask_test
#SBATCH --time=72:00:00
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --mail-type=end
#SBATCH --mail-user=szhan121@jhu.edu

module load r/4.3.0
module load anaconda
conda activate model_calibration
#pip install -U pip setuptools
#pip install -r requirement.txt

# Run the calibration
python model_calibration_spqsp.py