#!/bin/bash -l
#$ -N k200_clustering
#$ -j y
#$ -o k200_clustering.qlog
#$ -pe omp 8
#$ -P camplab
#$ -l mem_per_core=8G
#$ -l h_rt=12:00:00

module load R/4.5.2

Rscript h5_clustering_k20.R