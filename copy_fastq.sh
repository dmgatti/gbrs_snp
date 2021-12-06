#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 1 # number of cores
#SBATCH --mem 8G # memory pool for all cores
#SBATCH --time 0-1:00 # time (D-HH:MM)

################################################################################
# Copy some Attie Islet FASTQ files from source directory to /fastscratch.
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-05
################################################################################

##### VARIABLES #####

# Source directory.
SRC_DIR=/projects/churchill-lab/data/Attie/DO_islet_RNA/fastq/15NGS-003-gac_CHU11309_CHU11557

# Sample directories.
SAMPLE_DIRS=(`ls ${SRC_DIR}`)

# Destination directory.
DEST_DIR=/fastscratch/dgatti/fastq

# Number of samples to copy.
NUM_SAMPLES=20

mkdir -p ${DEST_DIR}

for I in $(seq 0 $((${NUM_SAMPLES}-1)))
do
  SAMPLE=${SAMPLE_DIRS[${I}]}

  echo Processing ${SAMPLE}

  SAMPLE_DIR=${SRC_DIR}/${SAMPLE}/fastq

  FASTQ_FILES=(`ls ${SAMPLE_DIR}`)

  cat ${SAMPLE_DIR}/${FASTQ_FILES[0]} ${SAMPLE_DIR}/${FASTQ_FILES[1]} > ${DEST_DIR}/${SAMPLE}.fastq.gz

done
