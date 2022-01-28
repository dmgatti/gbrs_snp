#!/bin/bash
#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 1 # number of cores
#SBATCH --mem 8G # memory pool for all cores
#SBATCH --time 0-2:00 # time (D-HH:MM)

################################################################################
# Create QUILT reference files.
#
# Get Sanger variants that fall withing annotated transcripts. Use and EnsemblDB 
# to get the transcript positions and query the Sanger Mouse Genomes SNP/INDEL file.
# Output to sorted and indexed VCF.
# Convert VCF to IMPUTE hap/legend/sample format.
# Physical map
# Genetic map
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2022-01-16
################################################################################

##### VARIABLES #####

# Sanger file directory.
SANGER_DIR=/projects/omics_share/mouse/GRCm38/genome/variants/snps_indels/rel_2004_v7/

# Path to sanger SNP/INDEL file.
SANGER_FILE=${SANGER_DIR}/mgp_REL2005_snps_indels.vcf.gz

# Ensembl version (we'll get teh EnsemblDB from AnnotationHub).
ENSEMBL_VERSION=102

# Reference file directory.
REF_DIR=/projects/compsci/USERS/dgatti/data/gbrs_snp

# R container with Bioconductor tools.
CONTAINER=~/containers/gbrs_snp_r.sif

# Samtools container.
SAMTOOLS=~/containers/samtools_1.10.sif

# mmconvert container.
MMCONVERT=~/containers/mmconvert.sif

# Script directory.
SCRIPT_DIR=/projects/compsci/USERS/dgatti/projects/gbrs_snp/scripts/gbrs_snp

# VCF script to run.
VCF_SCRIPT=${SCRIPT_DIR}/create_reference_files.R

# Genetic map script using mmconvert.
GMAP_SCRIPT=${SCRIPT_DIR}/create_gmap_files.R

##### PROGRAM #####

module load singularity

# Create VCF and physical map files.
singularity exec ${CONTAINER} Rscript ${VCF_SCRIPT} ${SANGER_FILE} ${ENSEMBL_VERSION} ${REF_DIR}

# Convert chromosome VCFs to IMPUTE haplotype and legend format.

VCF_FILES=(`ls ${REF_DIR}/*.vcf.bgz`)

for F in ${VCF_FILES[@]}
do
  # ${F%%.*z} deletes everything from "." to "z" at the end of the string.
  singularity exec ${SAMTOOLS} bcftools convert --haplegendsample ${F%%.*z} ${F}
done

# Create genetic map file.
singularity exec ${MMCONVERT} Rscript ${GMAP_SCRIPT} ${REF_DIR}

