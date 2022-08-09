#!/bin/bash
#SBATCH --job-name=net
#SBATCH --nodes=1 --ntasks-per-node=4
#SBATCH --time=100:00:00   # HH/MM/SS
#SBATCH --output=net_output
#SBATCH --mem=64G

module purge
module load qiime/intel/1.9.1

make_otu_network.py \
-i ../biom_table_no_chloroplast/feature-table-no-chloroplast.biom \
-m ../metadata.txt \
-o otu_network
