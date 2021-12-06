################################################################################
# Pileup the BAM files and intersect them with the Sanger variants that occur
# in transcripts.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-05
################################################################################

##### LIBRARIES #####

library(Rsamtools)

source('check_file_access.R')
source('read_sanger_vars.R')

##### VARIABLES #####

args = commandArgs(trailing = TRUE)

if(length(args) != 3) {

  stop(paste('ERROR: Three arguments required: bam_file, sanger_file, out_prefix'))

} # if(length(args) != 3)

bam_file    = args[1]
sanger_file = args[2]
out_prefix  = args[3]

# Test code
#bam_file   = '/fastscratch/dgatti/rsem/Sample_100-GES15-05705.STAR.genome.bam_sorted'
#sanger_file = '/projects/compsci/USERS/dgatti/data/gbrs_snp/sanger_transcript_snps_indels_ens102_b38.tsv.bgz'
#out_prefix  = '/fastscratch/dgatti/counts/Sample_100-GES15-05705'

check_file_access(bam_file,            4)
check_file_access(sanger_file,         4)
check_file_access(dirname(out_prefix), 2)

##### PROGRAM #####

# Open the Sanger file.
sanger = TabixFile(sanger_file)

header = headerTabix(sanger)

# Get the chromosomes.
chr_names = header$seqnames
chr_names = chr_names[order(as.numeric(chr_names))]

# Get column names.
cn = strsplit(sub('^#', '', header$header[length(header$header)]), '\t')[[1]]

# Process one chromosome at a time.
for(chr in chr_names) {

  print(paste('Chromosome', chr))
  
  # Read in Sanger variants on this chromosome.
  gr = GRanges(seqnames = chr, IRanges(start = 0, end = 200e6))
  vars = read_sanger_vars(sanger, gr)
  
  # Retain only SNPs for now.
  vars = subset(vars, indel == FALSE)
  
  # Retain founder SNPs that are homozygous.
  homo_alleles = c('A/A', 'T/T', 'G/G', 'C/C', '-/-')
  geno_columns = c('A_J', 'C57BL_6J', '129S1_SvImJ', 'NOD_ShiLtJ', 
                   'NZO_HlLtJ', 'CAST_EiJ', 'PWK_PhJ', 'WSB_EiJ')
  prop_homo = rowMeans(matrix(as.matrix(vars[,geno_columns]) %in% homo_alleles,
                       nrow = nrow(vars)))
  vars = subset(vars, prop_homo == 1.0)
  
  print(paste('Found', nrow(vars), 'variants on chr', chr))
  
  # Pileup BAM file on this chromosome.
  # TBD: Switch over to Rsamtools::PileupFiles?
  sb_param = ScanBamParam(what = 'seq', which = gr)
  pu_param = PileupParam()
  
  pu = pileup(bam_file, scanBamParam = sb_param, pileupParam = pu_param)
  pu = GRanges(pu$seqnames, IRanges(pu$pos, width = 1), strand = pu$strand,
               mcols = pu[,c('nucleotide', 'count')])
  colnames(mcols(pu)) = sub('^mcols\\.', '', colnames(mcols(pu)))
  
  vars = GRanges(vars$chr, IRanges(vars$pos, width = 1),
                 mcols = vars[,c('marker', 'ref', 'alt', 'tx_id', 'tx_strand', 'gene_id',
                                'indel', 'A_J', 'C57BL_6J', '129S1_SvImJ', 
                                'NOD_ShiLtJ', 'NZO_HlLtJ', 'CAST_EiJ', 'PWK_PhJ',
                                'WSB_EiJ')])
  colnames(mcols(vars)) = sub('^mcols\\.', '', colnames(mcols(vars)))

  ol = findOverlaps(pu, vars)
  
  pu = as.data.frame(pu[queryHits(ol)])
  pu = pu[,c('seqnames', 'start', 'strand', 'nucleotide', 'count')]
  vars = as.data.frame(vars[subjectHits(ol)])
  colnames(vars) = sub('^X', '', colnames(vars))
  vars = vars[,c('marker', 'ref', 'alt', 'tx_id', 'tx_strand', 'gene_id',
                 'indel', 'A_J', 'C57BL_6J', '129S1_SvImJ', 'NOD_ShiLtJ', 
                 'NZO_HlLtJ', 'CAST_EiJ', 'PWK_PhJ', 'WSB_EiJ')]
  counts = cbind(pu, vars)
    
  rm(pu, vars)
  gc()
  
  print(paste('Writing counts for', basename(out_prefix)))  
  
  counts = counts[,c('marker', 'seqnames', 'start', 'ref', 'alt', 'strand', 
                 'nucleotide', 'count', 'tx_id', 'tx_strand', 'gene_id',
                 'indel', 'A_J', 'C57BL_6J', '129S1_SvImJ', 'NOD_ShiLtJ', 
                 'NZO_HlLtJ', 'CAST_EiJ', 'PWK_PhJ', 'WSB_EiJ')]
  colnames(counts)[2:3] = c('chr', 'pos')
  
  out_file = paste0(out_prefix, '_counts.csv')
  write.table(counts, file = out_file, quote = FALSE, sep = ',', row.names = FALSE,
              col.names = chr == '1', append = chr != '1')
  
} # for(chr)

