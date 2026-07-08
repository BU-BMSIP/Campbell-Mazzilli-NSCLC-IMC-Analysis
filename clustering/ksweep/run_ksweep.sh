#!/bin/bash -l
#$ -N ksweep_broad_clustering
#$ -j y
#$ -o ksweep_broad_clustering.qlog
#$ -pe omp 8
#$ -P camplab
#$ -l mem_per_core=8G
#$ -l h_rt=12:00:00

module load R

Rscript ksweep.R
