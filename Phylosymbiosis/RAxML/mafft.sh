#!/bin/bash

#SBATCH --job-name=mafft_ASV
#SBATCH --output=mafft_ASV.log
#SBATCH --error=mafft_ASV.err

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G

#SBATCH --account=gwdef
#SBATCH --qos=high_prio

#SBATCH --partition=krypton

mafft --retree 2 --maxiterate 2 --thread 16 ASVs_clean_v138.2.fa > ASVs_clean_v138.2_aligned.fa
