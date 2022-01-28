#!/bin/bash
#SBATCH --qos batch
#SBATCH --partition compute
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 10 # number of cores
#SBATCH --mem 64G # memory pool for all cores
#SBATCH --time 0-8:00 # time (D-HH:MM)
#SBATCH --array 01-02

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

# Directory containing Sanger and reference files.
SANGER_DIR=/projects/compsci/USERS/dgatti/data/gbrs_snp

# Path to Sanger transcript SNP file for current chromosome.
SANGER_FILE=${SANGER_DIR}/sanger_chr${CHR}.vcf.bgz

# Reference haplotype file for current chromosome.
REF_HAP_FILE=${SANGER_DIR}/sanger_chr${CHR}.hap.gz

# Reference legend file.
REF_LEGEND_FILE=${SANGER_DIR}/sanger_chr${CHR}.legend.gz

# SNP physical position file for QUILT.
SNP_POS_FILE=${BASE_DIR}/snp_pos_chr${CHR}.txt

# SNP genetic map file for QUILT.
GEN_MAP_FILE=${SANGER_DIR}/sanger_chr${CHR}_gmap.txt.gz

# Number of cores to use. Must match --ntasks in SBATCH header above.
NUM_CORES=10

# Number of generations of outbreeding.
NUM_GEN=30

# Output directory.
OUT_DIR=${BASE_DIR}/chr${CHR}

# Output VCF file name.
VCF_FILE=${OUT_DIR}/quilt_chr${CHR}.vcf.gz

# QUILT container. This should have QUILT.R in the $PATH.
CONTAINER=~/containers/quilt.sif


##### PROGRAM #####

# Create output directory.
mkdir -p ${OUT_DIR}

module load singularity

# Write out BAM file list.
echo `ls -1 ${BAM_DIR}/*.STAR.genome.bam_sorted` > ${BAM_LIST_FILE} 

# Write out SNP position file.
# -P in grep to use perl syntax and catch the \t.
zcat ${SANGER_FILE} | sed 1,50d | cut -f 1,2,4,5 | grep -P "^${CHR}\t" > ${SNP_POS_FILE}

singularity exec ${CONTAINER} QUILT.R  \
                              --chr=${CHR} \
                              --regionStart=5000000 \
                              --regionEnd=7000000 \
                              --buffer=500000 \
                              --bamlist=${BAM_LIST_FILE} \
                              --posfile=${SNP_POS_FILE} \
                              --reference_haplotype_file=${REF_HAP_FILE} \
                              --reference_legend_file=${REF_LEGEND_FILE} \
                              --genetic_map_file=${GEN_MAP_FILE} \
                              --output_gt_phased_genotypes=TRUE \
                              --outputdir=${OUT_DIR} \
                              --output_filename=${VCF_FILE} \
                              --nGen=${NUM_GEN} \
                              --nCores=${NUM_CORES}

#                              --genfile=gen.txt \
