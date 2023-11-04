#!/bin/bash
#SBATCH --job-name Make_mask
#SBATCH --partition=kicp
#SBATCH --account=kicp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=1:00:00
#SBATCH --output=/home/dhayaa/DECADE/foreground_mask/make_mask.log
#SBATCH --mail-user=dhayaa@uchicago.edu
#SBATCH --mail-type=BEGIN,END


cd /home/dhayaa/DECADE/foreground_mask


#0.09655683845281593
for ext in 0.4 0.3 0.2 0.1 0.05
do
    out=/project/chihway/dhayaa/DECADE/Foreground_Masks/GOLD_Ext${ext}_Star5.fits
    python -u Runner/AllSkyRunner.py --max_extinction ${ext} --max_stardensity 5 --output_file_path ${out}
done


#3.044570751190183
for star in 2 3 5 10 20
do
    out=/project/chihway/dhayaa/DECADE/Foreground_Masks/GOLD_Ext0.2_Star${star}.fits
    python -u Runner/AllSkyRunner.py --max_extinction 0.2 --max_stardensity ${star} --output_file_path ${out}
done
