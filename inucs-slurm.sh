#!/bin/bash -l
#SBATCH --job-name=inucs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=100
#SBATCH --ntasks=100
#SBATCH --mem=600G
#SBATCH --output=slurm-4DNFI1O6IL1Q.pairs.log-20210909

echo "Server $HOSTNAME. Starting"
echo 

conda activate inucs

time python3.9 ./inucs.py prepare --multiprocessing 100 data/H1_chromosome.txt data/H1_nucleosome.txt data/4DNFI1O6IL1Q.pairs -d data/4DNFI1O6IL1Q.pairs.inucs

time python3.9 ./inucs.py plot data/4DNFI1O6IL1Q.pairs.inucs chr5 62380000 62440000 --save

echo 
echo "Server $HOSTNAME. Finished"
