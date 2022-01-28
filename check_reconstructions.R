################################################################################
# Examine the GBRS-SNP qtl2 haplotype reconstructions.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-06
################################################################################

##### LIBRARIES #####

library(tidyverse)
library(qtl2convert)
library(qtl2)

##### FUNCTIONS #####

# Return the maximum probability for each sample at each marker.
maxprob = function(genoprobs) {

  lapply(genoprobs, apply, c(1, 3), max)

} # maxprob()

# Get the number of crossovers.
getxo = function(geno) {

  tmp = lapply(geno, function(z) {
                            z[,-ncol(z)] - z[,-1]
                          })
  tmp = sapply(tmp, function(z) {
                      rowSums(z != 0) + 1
                    })
  rowSums(tmp)

} # getxo()

##### VARIABLES #####

min_coverage = '200'

geno_dir    = paste0('/fastscratch/dgatti/genotypes_', min_coverage)
results_dir = '/projects/compsci/USERS/dgatti/projects/gbrs_snp/results'
in_dir      = file.path(results_dir, 'qtl2')
#out_dir     = file.path(results_dir, 'genotype_intersect')
out_dir     = file.path(results_dir, 'genotype_union')

pmap_file = file.path(in_dir, 'pmap.csv')
gp_file = file.path(in_dir, paste0('genoprobs_36state_', min_coverage,'.rds'))
ap_file = file.path(in_dir, paste0('genoprobs_8state_', min_coverage,'.rds'))

samples = sub('_genotypes.csv$', '', dir(geno_dir))

##### PROGRAM #####

# Get number of variants per sample.
files = dir(geno_dir, pattern = '_genotypes.csv', full.names = TRUE)
num_geno = lapply(files, scan, what = 'char', sep = '\n')
num_geno = setNames(sapply(num_geno, length), samples)
names(num_geno) = gsub('-', '.', names(num_geno))

# Read in 36 state genoprobs.
genoprobs = readRDS(gp_file)

# Setting minprob low to get the maximum call.
gt = maxmarg(genoprobs, minprob = 0.01, quiet = FALSE)

# Get the maximum probability for each sample & marker.
max_prob = maxprob(genoprobs)
mean_prob  = sapply(max_prob, apply, 1, mean)
prob_gt_95 = sapply(max_prob, apply, 1, function(z) { mean(z > 0.95) })

# Get the number of crossovers.
xo = count_xo(gt, quiet = FALSE)
xo = rowSums(xo[,-ncol(xo)])

# Get the number of markers from the qtl2 pmap file.
num_markers = rowSums(sapply(genoprobs, dim))[3]

# Write out summary files.
num_geno = num_geno[names(xo)]
out_file = file.path(out_dir, paste0('hap_reconst_summary_min_cvrg_', min_coverage, '.csv'))
write.csv(data.frame(id     = names(xo), 
                     num_markers = num_markers,
                     num_var = num_geno,
                     num_xo = xo,
                     mean_probs = rowMeans(mean_prob),
                     mean_prob_gt_95 = rowMeans(prob_gt_95)),
          file = out_file, quote = FALSE, row.names = FALSE)

out_file = file.path(out_dir, paste0('hap_reconst_geno_min_cvrg_', min_coverage, '.rds'))
saveRDS(gt, out_file)
