################################################################################
# Given the count files from the BAM pileups, call genotypes at each variant
# and write out the marker, chr, position, and genotype.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-05
################################################################################

##### LIBRARIES #####

library(Rsamtools)

source('check_file_access.R')

##### VARIABLES #####

args = commandArgs(trailing = TRUE)

if(length(args) != 3) {

  stop(paste('ERROR: Three arguments required: in_file, min_covarage, out_path'))

} # if(length(args) != 3)

in_file      = args[1]
min_coverage = as.numeric(args[2])
out_path     = args[3]

# Test code
#in_file       = '/fastscratch/dgatti/counts/Sample_100-GES15-05705_counts.csv'
#min_coverage = 10
#out_path     = '/fastscratch/dgatti/genotypes'

check_file_access(in_file,           4)
check_file_access(sanger_file,       4)
check_file_access(dirname(out_path), 2)

if(is.na(min_coverage)) {

  stop(paste('ERROR: Minimum coverage  must be an integer. Found:', args[2]))

} # if(is.na(min_coverage))

##### PROGRAM #####

print(basename(in_file))

# Read in the counts file.
counts = read.csv(in_file)

print(paste('Read in', nrow(counts), 'rows.'))

# Change the count column to numeric. There are sometimes '-' calls, which 
# produces letters.
counts$count = as.numeric(counts$count)
counts       = subset(counts, !is.na(counts$count))

# Filter by minimum coverage.
counts = subset(counts, count >= min_coverage)
stopifnot(min(counts$count, na.rm = TRUE) == min_coverage)

print(paste(nrow(counts), 'after filtering by minimum coverage', min_coverage))

# Filter to retain genotypes with canonical allele calls that match the Sanger alleles.
# This also filters out rows for which there are more than one alt allele.
counts = subset(counts, nucleotide == ref | nucleotide == alt)

print(paste(nrow(counts), 'after filtering for canonical alleles.'))

# TBD: Handle females and males differently becuase of differences in number
#      of copies.

# NOTE: I would have expected most of the pileup strands to match the transcript
# strans. They don't. Maybe this is becuase we're sequencing cDNA?
#> table(counts$strand)
#     -      +
#200941 197029
#> table(counts$tx_strand)
#     -    -;+      +    +;-
#180982   8671 183836  24481
#> mean(counts$strand == counts$tx_strand)
#[1] 0.0206699

# We still have the transcript, gene, and founder alleles if we need them.
# For now, remove them.
counts = counts[,c('marker', 'chr', 'pos', 'strand', 'nucleotide', 'count')]

# Call the genotypes.
# TBD: In the long run, we should write a custom HMM to weight the genotype
# probabilities by the read depth.
counts = split(counts, counts$marker)
# Retain four columns.
counts = lapply(counts, function(z) {
                          z[, c('marker', 'chr', 'pos', 'nucleotide')]
                        })
# Get genotype. 
counts = lapply(counts, function(z) {
                          z$nucleotide[1] = paste(sort(unique(z$nucleotide)), collapse = '')
                          z
                        })
# Keep the first row of each element.
counts = sapply(counts, function(z) {
                          z[1,]
                        })
# Adjust data format. (There may be an easier way)
counts = t(counts)
rownames(counts) = NULL
counts = data.frame(counts)
counts$marker     = unlist(counts$marker)
counts$chr        = unlist(counts$chr)
counts$pos        = unlist(counts$pos)
counts$nucleotide = unlist(counts$nucleotide)
colnames(counts)[4] = 'gt'

# Make single genotype values homozygous.
# TBD: Handle Chr X correctly.
wh = which(nchar(counts$gt) == 1)
counts$gt[wh] = paste0(counts$gt[wh], counts$gt[wh])

stopifnot(all(nchar(counts$gt) <= 2))

# Write out the file.
out_file = file.path(out_path, sub('counts', 'genotypes', basename(in_file)))
print(paste('Writing', out_file))

write.table(counts, file = out_file, sep = ',', quote = FALSE, row.names = FALSE)



# Plot counts vs het/homo.
# Hets have higher counts.
#tmp = split(counts, counts$marker)
#tmp = sapply(tmp, function(z) {
#                         gt = paste0(unique(z$nucleotide), collapse = '')
#                         c(gt, sum(z$count))
#                       })
#tmp = t(tmp)
#wh1 = which(nchar(tmp[,1]) == 1)
#wh2 = which(nchar(tmp[,1]) > 1)
#tmp[wh1,1] = 'homo'
#tmp[wh2,1] = 'het'
