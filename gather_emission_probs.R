################################################################################
# Given the count files from the BAM pileups, estimate the emission 
# probabilities from each SNP.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-09
################################################################################

##### LIBRARIES #####

library(Rsamtools)

source('check_file_access.R')

##### VARIABLES #####

args = commandArgs(trailing = TRUE)

if(length(args) != 3) {

  stop(paste('ERROR: Three arguments required: in_path, min_covarage, out_path'))

} # if(length(args) != 3)

in_path      = args[1]
min_coverage = as.numeric(args[2])
out_path     = args[3]

# Test code
#in_path       = '/fastscratch/dgatti/counts'
#min_coverage = 1
#out_path     = '/fastscratch/dgatti/emission_probs'

check_file_access(in_file,           4)
check_file_access(dirname(out_path), 2)

if(is.na(min_coverage)) {

  stop(paste('ERROR: Minimum coverage  must be an integer. Found:', args[2]))

} # if(is.na(min_coverage))

##### PROGRAM #####

# Get the counts files.
count_files = dir(in_path, pattern = '_counts\\.csv$')

# Marker positions and alleles. 
markers = NULL

for(f in count_files) {

  filename = file.path(in_path, f)
  sample = sub('_counts\\.csv', '', f)
  
  # Read counts file.
  counts = read.csv(filename)
  # Filter to retain minimum coverage.
  counts = subset(counts, count >= min_coverage)

  
} # for(f)
