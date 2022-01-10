#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 1 # number of cores
#SBATCH --mem 8G # memory pool for all cores
#SBATCH --time 0-1:00 # time (D-HH:MM)
#SBATCH --array=000-128

################################################################################
# Pileup counts in genome BAM files (sorted & indexed) at the Sanger transcript
# variants. Output a file with counts and loci.
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-05
################################################################################

##### VARIABLES #####

# Full path to directory containing BAM files.
BAM_DIR=/fastscratch/dgatti/rsem

# Full path of BAM files in BAM_DIR.
BAM_FILES=(`ls ${BAM_DIR}/*STAR.genome.bam_sorted`)

# Full path to current BAM file.
BAM_FILE=${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}

# Path to Sanger transcript SNP/INDEL file.
SANGER_FILE=/projects/compsci/USERS/dgatti/data/gbrs_snp/sanger_transcript_snps_indels_ens102_b38.tsv.bgz

# Sample ID for this sample.
SAMPLE=${BAM_FILE/${BAM_DIR}\//}
SAMPLE=${SAMPLE/.STAR.genome.bam_sorted/}

# Full path to output directory.
OUT_DIR=/fastscratch/dgatti/counts

# Output file prefix for this sample.
OUT_FILE=${OUT_DIR}/${SAMPLE}

# R Bioconductor container.
CONTAINER=~/containers/gbrs_snp_r.sif

# R script to run.
RSCRIPT=get_variant_counts.R

##### PROGRAM #####

mkdir -p ${OUT_DIR}

touch ${BAM_DIR}/*

module load singularity

echo Processing ${SAMPLE}

singularity exec ${CONTAINER} Rscript ${RSCRIPT} ${BAM_FILE} ${SANGER_FILE} ${OUT_FILE}

