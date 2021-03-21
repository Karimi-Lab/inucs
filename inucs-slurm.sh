#!/bin/bash -l
#SBATCH --job-name=inucs
#SBATCH --mem=44G
#SBATCH --output=slurm-4DNFI1O6IL1Q.pairs.gz.log
# From man sbatch: SBATCH --mem=0
# "A memory size specification of zero is treated as a special case and grants the job access to all of the memory on each node."
echo "Server $HOSTNAME. Starting"
echo 

conda activate inucs

# ./inucs.py prepare --help
time python3.8 ./inucs.py prepare data/H1_chromosome.txt data/H1_nucleosome.txt data/4DNFI1O6IL1Q.pairs.gz -d data/4DNFI1O6IL1Q.pairs.gz.inucs

time python3.8 ./inucs.py plot data/4DNFI1O6IL1Q.pairs.gz.inucs chr2 27325465 27442377 --save

echo 
echo "Server $HOSTNAME. Finished"
