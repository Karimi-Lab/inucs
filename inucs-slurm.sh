#!/bin/bash -l
#SBATCH --job-name=inucs
#SBATCH --mem-per-cpu=5G
# From man sbatch: # SBATCH --mem=0
# "A memory size specification of zero is treated as a special case and grants the job access to all of the memory on each node."
echo "Server $HOSTNAME. Starting"
echo 

conda activate inucs

#./inucs.py prepare --help
time ./inucs.py prepare data/H1_chromosome.txt data/H1_nucleosome.txt data/4DNFI1O6IL1Q.pairs.gz -d data/4DNFI1O6IL1Q.pairs.gz.inucs

echo 
echo "Server $HOSTNAME. Finished"
