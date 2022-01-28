#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 10 # number of cores
#SBATCH --mem 64G # memory pool for all cores
#SBATCH --time 0-8:00 # time (D-HH:MM)
#SBATCH --array 01-20

################################################################################
# Run QUILT on the BAM files from the RNAseq alignment.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-01-11
################################################################################

##### VARIABLES #####

# Base directory where QUILT configuration files are located.
BASE_DIR=/fastscratch/dgatti/quilt

# Directory where BAM files are located..
BAM_DIR=/fastscratch/dgatti/rsem

# File containing list of BAM files.
BAM_LIST_FILE=${BASE_DIR}/bam_list.txt

# Chromosome to run.
CHR=${SLURM_ARRAY_TASK_ID}

if [ ${CHR} -eq 20 ]
then
  CHR="X"
fi

# Path to Sanger transcript SNP/INDEL file.
SANGER_FILE=/projects/compsci/USERS/dgatti/data/gbrs_snp/sanger_transcript_snps_indels_ens102_b38.tsv.bgz

# SNP position file for STITCH input.
SNP_POS_FILE=${BASE_DIR}/snp_pos_chr${CHR}.txt

# Number of founder strains.
NUM_FOUNDERS=8

# Number of cores to use. Must match --ntasks in SBATCH header above.
NUM_CORES=10

# Number of generations of outbreeding.
NUM_GEN=30

# Output directory.
OUT_DIR=${BASE_DIR}/chr${CHR}

# Output VCF file name.
VCF_FILE=${OUT_DIR}/stitch_chr${CHR}.vcf.gz

# QUILT container. This should have QUILT.R in the $PATH.
CONTAINER=~/containers/quilt.sif


##### PROGRAM #####

# Create output directory.
mkdir -p ${OUT_DIR}

module load singularity

# Write out BAM file list.
cat `ls -1 ${BAM_DIR}/*.STAR.genome.bam_sorted` > ${BAM_LIST_FILE} 

# Write out SNP position file.
# -P in grep to use perl syntax and catch the \t.
zcat ${SANGER_FILE} | sed 1,8d | cut -f 2,3,4,5 | grep -P "^${CHR}\t" > ${SNP_POS_FILE}

singularity exec ${CONTAINER} QUILT.R  \
                              --chr=${CHR} \
                              --regionStart=5000000 \
                              --regionEnd=6000000 \
                              --buffer=10000 \
                              --bamlist=${BAM_LIST_FILE} \
                              --posfile=${SNP_POS_FILE} \
                              --output_haplotype_dosages=TRUE \
                              --outputdir=${OUT_DIR} \
                              --output_filename=${VCF_FILE} \
                              --K=${NUM_FOUNDERS} \
                              --nGen=${NUM_GEN} \
                              --nCores=${NUM_CORES}

#                              --genfile=gen.txt \
# --genetic_map_file
# --reference_haplotype_file

