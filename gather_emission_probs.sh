#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 1 # number of cores
#SBATCH --mem 64G # memory pool for all cores
#SBATCH --time 0-2:00 # time (D-HH:MM)

################################################################################
# Given the piled up counts files from each sample, gather the genotypes across
# all samples and produce emission probabilities for each marker.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-10
################################################################################

##### VARIABLES #####

# Full path to input file directory.
IN_DIR=/fastscratch/dgatti/counts

# Listing of all input files.
IN_FILES=(`ls ${IN_DIR}/*_counts.csv`)

# Minimum read depth to retain a variant call.
# TBD: Add some mapping quality metric?
MIN_COVERAGE=1

# Full path to output directory.
OUT_DIR=/fastscratch/dgatti/emission_probs

# R Bioconductor container.
CONTAINER=~/containers/gbrs_snp_r.sif

# R script to run.
RSCRIPT=gather_emission_probs.R

##### PROGRAM #####

mkdir -p ${OUT_DIR}

module load singularity

singularity exec ${CONTAINER} Rscript ${RSCRIPT} ${IN_FILE} ${MIN_COVERAGE} ${OUT_DIR}



