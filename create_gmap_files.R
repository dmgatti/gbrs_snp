################################################################################
# Create genetic map files, one per chromosome, for QUILT.
# Assumes that gzipped IMPUTE legend files exist in ref_dir.
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2022-01-21
################################################################################

##### LIBRAIRES #####

library(mmconvert)
source('check_file_access.R')

##### VARIABLES #####

args = commandArgs(trailing = TRUE)

if(length(args) != 1) {

  stop(paste('ERROR: One argument required: ref_dir'))

} # if(length(args) != 1)

ref_dir = args[1]

# ref_dir = '/projects/compsci/USERS/dgatti/data/gbrs_snp'

check_file_access(ref_dir, 2)

pmap_files = dir(ref_dir, pattern = 'legend.gz', full.names = TRUE)

##### PROGRAM #####

chromosomes = gsub('^sanger_chr|\\.legend\\.gz$', '', basename(pmap_files))

for(i in seq_along(pmap_files)) {
  
  chr = chromosomes[i]

  pmap = read.delim(pmap_files[i], header = TRUE, sep = ' ')
  pmap = data.frame(chr    = chr,
                    pos    = pmap$position,
                    marker = pmap$id)
  gmap = mmconvert(positions = pmap, input_type = 'bp')
  
  if(chr != 'X') {

    gmap = data.frame(pos  = pmap$pos,
                      rate = 0,
                      cM   = gmap$cM_coxV3_ave)

  } else {

    gmap = data.frame(pos  = pmap$pos,
                      rate = 0,
                      cM   = gmap$cM_coxV3_female)
  
  } # else
  
  # g[iRow, 3]  = g[iRow - 1, 3] + ( (g[iRow, 1] - g[iRow - 1, 1]) * g[iRow - 1, 2])
  gmap$cM = gmap$cM - min(gmap$cM)
  gmap$rate[1:(nrow(gmap) - 1)] = diff(gmap$cM) / diff(gmap$pos * 1e-6)
  colnames(gmap) = c('position', 'COMBINED_rate.cM.Mb.', 'Genetic_Map.cM.')

  # Temporary filtering.
  gmap = subset(gmap, position <= 8e6)

  out_file = file.path(ref_dir, paste0('sanger_chr', chr, '_gmap.txt'))
  write.table(gmap, file = out_file, sep = ' ', row.names = FALSE, 
              col.names = FALSE, quote = FALSE)

  system(paste('gzip -f', out_file))

} # for(i)


