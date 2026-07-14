#!/bin/bash -l
#$ -N k200_clustering
#$ -j y
#$ -o k200_clustering.qlog
#$ -pe omp 8
#$ -P camplab
#$ -l mem_per_core=8G
#$ -l h_rt=12:00:00

module load R

Rscript /restricted/projectnb/camplab/projects/20220504_Suzuki_LN/TumorIMC/analysis/clustering/h4_k200_clusters/k200.R
