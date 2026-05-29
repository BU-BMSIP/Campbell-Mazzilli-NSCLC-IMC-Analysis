#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -P camplab
#$ -N runSteinbock
#$ -o qsub_steinbock.log


echo "=========================================================="
echo "Starting on       : $(date)"
echo "Running on node   : $(hostname)"
echo "Current job ID    : $JOB_ID"
echo "Current job name  : $JOB_NAME"
echo "Task index number : $TASK_ID"
echo "=========================================================="


# load modules
module load steinbock
module load mpfr gcc R/4.2.1  gcc/12.2.0 R/4.2.1 zlib/1.2.13 sqlite3/3.37.2 geos/3.11.1 gdal/3.6.4 proj/9.2.0  python3/3.10.12 texlive pandoc

# run script for pixel filtration
#roiDirectory=""
#outputDirectory=""

Rscript scripts/runPixelFiltration.R ${roiDirectory} ${outputDirectory}

# convert .txt file in the 'raw' folder to TIFF and filter hot pixels
steinbock preprocess imc images --hpf 100

# deep learning-based cell segmentation
steinbock segment deepcell --minmax -o masks_deepcell

# measurement of the intensity for each cell and each marker
steinbock measure intensities --masks masks_deepcell
steinbock measure regionprops --masks masks_deepcell
steinbock measure neighbors --masks masks_deepcell --type expansion --dmax 4

# export the ome.tiff 
steinbock export ome

# export the single cell data in different formats
steinbock export histocat --masks masks_deepcell
steinbock export csv intensities regionprops -o cells.csv
steinbock export csv intensities regionprops --no-concat -o cells_csv
steinbock export fcs intensities regionprops -o cells.fcs
steinbock export fcs intensities regionprops --no-concat -o cells_fcs
steinbock export anndata --intensities intensities --data regionprops --neighbors neighbors -o cells.h5ad
steinbock export anndata --intensities intensities --data regionprops --neighbors neighbors --no-concat -o cells_h5ad
steinbock export graphs --data intensities

