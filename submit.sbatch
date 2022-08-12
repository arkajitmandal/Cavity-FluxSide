#!/bin/bash
#SBATCH -p action -A action
#SBATCH -J avemaria 
#SBATCH -o my_output_%j
#SBATCH --mem-per-cpu=2GB
#SBATCH -t 1:00:00
##SBATCH -N 1 
python3 main.py 10 10