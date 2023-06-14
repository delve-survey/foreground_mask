#!/bin/bash
#SBATCH --job-name Make_mask
##SBATCH --partition=broadwl
#SBATCH --partition=chihway
#SBATCH --account=pi-chihway
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=28
#SBATCH --ntasks-per-node=40
#SBATCH --time=12:00:00
#SBATCH --output=/home/dhayaa/Desktop/DECADE/foreground_mask/make_mask.log
#SBATCH --mail-user=dhayaa@uchicago.edu
#SBATCH --mail-type=BEGIN,END


cd /home/dhayaa/Desktop/DECADE/foreground_mask

python -u Runner/AllSkyRunner.py
