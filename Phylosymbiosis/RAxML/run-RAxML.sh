#!/bin/bash

#SBATCH --job-name=RAxML_ASV
#SBATCH --output=RAxML_ASV.log
#SBATCH --error=RAxML_ASV.err

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G

#SBATCH --account=gwdef
#SBATCH --qos=high_prio

#SBATCH --partition=krypton


~/RAxML/standard-RAxML-8.2.13/raxmlHPC-PTHREADS-AVX -f d -m GTRGAMMA -p 16647 -s ASVs_clean_v138.2_aligned.fa -n ASV_tree_v138.2 -T 32

#~/RAxML/standard-RAxML-8.2.13/raxmlHPC-PTHREADS-AVX  -f a -m GTRGAMMA -p 16647 -x 58421 -# 100 -s ASVs_clean_v138.2_aligned.fa -n ASV_tree_v138.2 -T 32



