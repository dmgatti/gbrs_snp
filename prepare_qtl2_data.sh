#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 1 # number of cores
#SBATCH --mem 32G # memory pool for all cores
#SBATCH --time 0-1:00 # time (D-HH:MM)

################################################################################
# Gather the sample & founder genotypes and write out their values at the 
# union of all markers. Also write out physical and genetic maps.
# I'm currently hard-coding the covariate & phenotype data.
# TBD: How to handle phenotype & covariate data.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-05
################################################################################

##### VARIABLES #####

# Full path to the genotype files.
IN_DIR=/fastscratch/dgatti/genotypes_200

# Full path to Sanger transcript variant file.
SANGER_FILE=/projects/compsci/USERS/dgatti/data/gbrs_snp/sanger_transcript_snps_indels_ens102_b38.tsv.bgz

# Full path to the output directory to contain the qtl2 files.
OUT_DIR=/projects/compsci/USERS/dgatti/projects/gbrs_snp/results/qtl2

# R Bioconductor container.
CONTAINER=~/containers/gbrs_snp_r.sif

# R script to run.
RSCRIPT=prepare_qtl2_data.R

##### PROGRAM #####

mkdir -p ${OUT_DIR}

module load singularity

singularity exec ${CONTAINER} Rscript ${RSCRIPT} ${IN_DIR} ${SANGER_FILE} ${OUT_DIR}

