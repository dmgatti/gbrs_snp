################################################################################
# Given the qtl2 files in the qtl2 directory, read in the cross, calculate
# genoprobs and allele probs and write them out to *.rds files.
# Arguments:
# json_file: Full path to the qtl2 JSON control file.
# num_cores: Integer that is the number of cores to use.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-05
################################################################################

##### LIBRARIES #####

library(qtl2)

source('check_file_access.R')

##### VARIABLES #####

args = commandArgs(trailing = TRUE)

if(length(args) != 2) {

  stop(paste('ERROR: Two arguments required: json_file, num_cores'))

} # if(length(args) != 1)

json_file = args[1]
num_cores = as.numeric(args[2])

# Test code
#json_file = '/projects/compsci/USERS/dgatti/projects/gbrs_snp/results/qtl2/control.json'

check_file_access(json_file, 4)
check_file_access(dirname(json_file), 2)

if(is.na(num_cores)) {

  stop(paste('ERROR" num_cores must be an integer.', num_cores))

} # if(is.na(num_cores))

# Output diretory.
# TBD: Allow user to write probs somewhere else.
out_dir = dirname(json_file)

##### PROGRAM #####

# Read in cross.
cross = read_cross2(json_file, quiet = FALSE)

print(summary(cross))

print('Calcuating genoprobs')

genoprobs = calc_genoprob(cross, error_prob = 0.01, quiet = FALSE, cores = num_cores)
saveRDS(genoprobs, file.path(out_dir, 'genoprobs_36state.rds'))

print('Calculating alleleprobs')

alleleprobs = genoprob_to_alleleprob(genoprobs, quiet = FALSE, cores = num_cores)
saveRDS(alleleprobs, file.path(out_dir, 'genoprobs_8state.rds'))


