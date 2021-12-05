# gbrs_snp
Genotyping by RNAseq using SNPs in Transcripts

This is currently a prototype. It needs to be converted to a pipeline for 
push-button operation.

# Steps in Analysis

Each script has a VARIABLES section at the top. The ones that need to be set are listed in the input section in each step.

1. Get Sanger SNPs & indels in DO founders that overlap with Ensembl transcripts. 
   Inputs:
   * SANGER_FILE: Full path to Sanger SNP/Indel VCF.
   * ENSEMBL_VERSION: Number of Ensembl version.
   * OUT_FILE: Full path to Sanger transcript variant output file.

   Output is a bgzipped, Tabix indexed file containing the SNPs & indels.
   
   Run create_sanger_variant_file.sh
     which calls create_sanger_variant_file.R

2. Create STAR/RSEM index.
   Inputs:
   * FASTA_FILE: Full path to mouse genome FASTA file.
   * GTF_FILE: Full path to Ensembl GTF file. Must match Ensembl version in step 1.
   * INDEX_DIR: Full path to where the STAR/RSEM index will be written.
   * READ_LENGTH: Read length of target FASTQ files that will align to this index.

   Output is a set of STAR/RSEM genome indices.
   
   Run rsem_index_genome.sh.

3. Align FASTQ files to transcriptome.
   Inputs:
   * INDEX_DIR: Full path to the STAR/RSEM index.
   * FASTQ_DIR: Full path to the FASTQ files.
   * OUT_DIR: Full path to directory where aligned files are written.
   * TEMP_DIR: Full path to a temporary working directory for STAR/RSEM.
   * DEST_DIR: Full path to the directory where the gene count files should be stored.

   Output is a set of STAR/RSEM transcriptome alignments. Sorted, indexed genome BAM files are produced, which are used in succcessive steps. Gene count files are produced, which could be used for gene quantification. Other DO-specific tools may be superior.

   Run rsem_quantify.sh.

4. Pileup BAM files and extract the counts at each variant.
   Inputs:
   * BAM_DIR: Full path to directory where sorted, indexed BAM files are.
   * SANGER_FILE: Full path to Sanger transcript variants file produced in step 1.
   * OUT_FILE: Full path to output directory, with sample prefix.

   The script will run each file in the SLURM ARRAY, so you'll need to adjust the array paramenter in the SBATCH header.

   Output is a set of *.csv files which contain the counts at each variant.

   Run get_variant_counts.sh
     which calls get_variant_counts.R

5. Select variants and call genotypes at each variants. We currently allow the user to filter by read depth. Perhaps some quality metric from the BAM file would be better or a good addition.
   Inputs:
   * IN_DIR: Full path to input file directory containing counts files from step 4.
   * MIN_COVERAGE: Integer that is the minimum coverage at which to call genotypes.
   * OUT_DIR: Full path to output file directory that will contain genotype files.

   Output is a set of *.csv files with contain genotype calls for each sample
   The script will run each file in the SLURM ARRAY, so you'll need to adjust the array paramenter in the SBATCH header.

   Run call_genotypes.sh
     which calls call_genotypes.R

6. Produce qtl2 input files. Gather the genotypes for all samples, get founder genotype calls for each variant, write out map files, covariates, phenotypes, and the json file for qtl2.
   Inputs:
   * IN_DIR: Full path to input file directory containing genotype files from step 5.
   * SANGER_FILE: Full path to Sanger transcript variants file produced in step 1.
   * OUT_DIR: Full path to the output directory which will hold all files required for qtl2 input.

  Output is the set of *.csv files required for qtl2, and the *.json control file.
   * sample_geno.csv: Sample genotypes, markers in rows, samples in columns.
   * founder_geno.csv: Founder genotypes, markers in rows, samples in columns.
   * pmap.csv: Physical map, markers in rows. Positions in Mb.
   * gmap.csv: Genetic map, markers in rows. Positions in cM.
   * covar.csv: Covariates, samples in rows, covariates in columns.
   * pheno.csv: Phenotypes, samples in rows, covariates in columns.
   * control.json: 

   Run prepare_qtl2_data.sh
     which calls prepare_qtl2_data.R
