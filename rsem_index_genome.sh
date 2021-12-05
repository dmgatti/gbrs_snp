#!/bin/bash
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 48G # memory pool for all cores
#SBATCH -t 0-2:00 # time (D-HH:MM)

################################################################################
# Build RSEM/STAR genome index for GRCm38.
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-05
################################################################################

### VARIABLES ###

# Directory containing FASTA file for mouse genome.
FASTA_DIR=/projects/omics_share/mouse/GRCm38/genome/sequence/fasta/ensembl

# Name of mouse genome FASTA file (unzipped).
FASTA_FILE=${FASTA_DIR}/Mus_musculus.GRCm38.dna.primary_assembly.fa

# Directory containing Ensembl GTF file.
GTF_DIR=/projects/omics_share/mouse/GRCm38/transcriptome/annotation/ensembl/gtf

# File path to Ensembl GTF file (unzipped).
GTF_FILE=${GTF_DIR}/Mus_musculus.GRCm38.102.gtf

# Output path for STAR/RSEM index with prefix.
INDEX_DIR=/fastscratch/dgatti/index/rsem

# Read length of FASTQ files.
READ_LENGTH=101

# STAR/RSEM singularity container.
CONTAINER=~/containers/star_rsem.sif


### PROGRAM ###

module load singularity

echo 'Building genome reference'

singularity exec ${CONTAINER} rsem-prepare-reference \
                              --num-threads 8 \
                              --star \
                              --star-sjdboverhang $(($READ_LENGTH - 1)) \
                              --gtf ${GTF_FILE} \
                              ${FASTA_FILE} \
                              ${INDEX_DIR}


