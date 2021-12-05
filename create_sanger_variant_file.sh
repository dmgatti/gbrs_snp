#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 1 # number of cores
#SBATCH --mem 32G # memory pool for all cores
#SBATCH --time 0-2:00 # time (D-HH:MM)

################################################################################
# Create a file containing the Sanger variants that fall withing annotated
# transcripts. Use and EnsemblDB to get the transcript positions and query
# the Sanger Mouse Genomes SNP/INDEL file.
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-05
################################################################################

##### VARIABLES #####

# Sanger file directory.
SANGER_DIR=/projects/omics_share/mouse/GRCm38/genome/variants/snps_indels/rel_2004_v7/

# Path to sanger SNP/INDEL file.
SANGER_FILE=${SANGER_DIR}/mgp_REL2005_snps_indels.vcf.gz

# Ensembl version (we'll get teh EnsemblDB from AnnotationHub).
ENSEMBL_VERSION=102

# Full path to the output file.
OUT_FILE=/projects/compsci/USERS/dgatti/data/gbrs_snp/sanger_transcript_snps_indels_ens102_b38.tsv

# R container with Bioconductor tools.
CONTAINER=~/containers/bioconductor.sif

# R script to run.
RSCRIPT=/projects/compsci/USERS/dgatti/projects/gbrs_snp/scripts/create_sanger_variant_file.R

##### PROGRAM #####

module load singularity

singularity exec ${CONTAINER} Rscript ${RSCRIPT} ${SANGER_FILE} ${ENSEMBL_VERSION} ${OUT_FILE}

