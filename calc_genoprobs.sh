#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 10 # number of cores
#SBATCH --mem 64G # memory pool for all cores
#SBATCH --time 0-8:00 # time (D-HH:MM)

################################################################################
# Given the qtl2 files in the qtl2 directory, read in the cross, calculate
# genoprobs and allele probs and write them out to *.rds files.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-06
################################################################################

##### VARIABLES #####

# Full path to the output directory to contain the qtl2 files.
QTL2_DIR=/projects/compsci/USERS/dgatti/projects/gbrs_snp/results/qtl2/control.json

# Number of cores to use. Must match --ntasks argument in slurm header above.
NUM_CORES=10

# R Bioconductor container.
CONTAINER=~/containers/gbrs_snp_r.sif

# R script to run.
RSCRIPT=calc_genoprobs.R

##### PROGRAM #####

module load singularity

singularity exec ${CONTAINER} Rscript ${RSCRIPT} ${QTL2_DIR} ${NUM_CORES}

