#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 1 # number of cores
#SBATCH --mem 8G # memory pool for all cores
#SBATCH --time 0-1:00 # time (D-HH:MM)
#SBATCH --array=00-19

################################################################################
# Once the counts have been piled up at the Sanger variants, filter them and
# produce genotype calls at each variant.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-05
################################################################################

##### VARIABLES #####

# Full path to input file directory.
IN_DIR=/fastscratch/dgatti/counts

# Listing of all input files.
IN_FILES=(`ls ${IN_DIR}/*_counts.csv`)

# Full path to input file for this run.
IN_FILE=${IN_FILES[${SLURM_ARRAY_TASK_ID}]}

# Full path to Sanger transcript variant file.
SANGER_FILE=/projects/compsci/USERS/dgatti/data/gbrs_snp/sanger_transcript_snps_indels_ens102_b38.tsv

# Minimum read depth to retain a variant call.
# TBD: Add some mapping quality metric?
MIN_COVERAGE=10

# Full path to output directory.
OUT_DIR=/fastscratch/dgatti/genotypes

# R Bioconductor container.
CONTAINER=~/containers/bioconductor.sif

# R script to run.
RSCRIPT=call_genotypes.R

##### PROGRAM #####

mkdir -p ${OUT_DIR}

module load singularity

echo Processing `basename ${IN_FILE}`

singularity exec ${CONTAINER} Rscript ${RSCRIPT} ${IN_FILE} ${SANGER_FILE} ${MIN_COVERAGE} ${OUT_DIR}


