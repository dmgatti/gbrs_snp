#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 8 # number of cores
#SBATCH --mem 32G # memory pool for all cores
#SBATCH --time 0-2:00 # time (D-HH:MM)
#SBATCH --array=00-129

################################################################################
# Use RSEM & STAR to estimate gene counts aligned to mouse genome.
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-05
################################################################################

### VARIABLES ###

# Directory containing STAR/RSEM genome index.
INDEX_DIR=/fastscratch/dgatti/index/rsem

# FASTQ file directory.
FASTQ_DIR=/fastscratch/dgatti/fastq

# Listing of all FASTQ files.
FASTQ_FILES=(`ls ${FASTQ_DIR}`)

# FASTQ file for this run.
FASTQ_FILE=${FASTQ_DIR}/${FASTQ_FILES[${SLURM_ARRAY_TASK_ID}]}

# Output directory for results.
OUT_DIR=/fastscratch/dgatti/rsem

# Temporary working directory.
TEMP_DIR=/fastscratch/dgatti/temp

# Final destination for results files.
DEST_DIR=/projects/compsci/USERS/dgatti/projects/gbrs_snp/results/rsem

# STAR/RSEM singularity container.
RSEM_CONTAINER=~/containers/star_rsem.sif

# Samtools container.
SAMTL_CONTAINER=~/containers/samtools_1.10.sif

### PROGRAM ###

mkdir -p ${OUT_DIR}

module load singularity

SAMPLE=`basename ${FASTQ_FILE}`
SAMPLE=${SAMPLE/.fastq.gz/}

echo Aligning ${SAMPLE}

# BAM file.
BAM=${OUT_DIR}/${SAMPLE}.STAR.genome.bam

mkdir -p ${TEMP_DIR}/${SAMPLE}

echo ${TEMP_DIR}/${SAMPLE}
echo ${FASTQ_FILE}
echo ${INDEX_DIR}
echo ${OUT_DIR}/${SAMPLE}

# NOTE: Ran into an error when trying to get RSEM to output a genome BAM.
# rsem-tbam2gbam /fastscratch/dgatti/index/rsem /fastscratch/dgatti/rsem/Sample_100-GES15-05705.transcript.bam /fastscratch/dgatti/rsem/Sample_100-GES15-05705.genome.bam
#rsem-tbam2gbam: BamConverter.h:131: void BamConverter::process(): Assertion `cqname != qname' failed.

singularity exec ${RSEM_CONTAINER} rsem-calculate-expression \
                                   --num-threads 8 \
                                   --star \
                                   --star-gzipped-read-file \
                                   --star-output-genome-bam \
                                   --temporary-folder ${TEMP_DIR}/${SAMPLE} \
                                   ${FASTQ_FILE} \
                                   ${INDEX_DIR} \
                                   ${OUT_DIR}/${SAMPLE}

# Sort BAM
singularity exec ${SAMTL_CONTAINER} samtools sort -@ 8 -o ${BAM}_sorted ${BAM}

# Index BAM
singularity exec ${SAMTL_CONTAINER} samtools index ${BAM}_sorted
rm ${BAM}

