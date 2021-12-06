################################################################################
# Read in the Sanger variants in the requested interval.
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-05
################################################################################

library(Rsamtools)

# Arguments:
# sanger: TabixFile that is the Sanger transcript variant file.
# param: GRanges object with interval to query.
read_sanger_vars = function(sanger, param) {

  header = headerTabix(sanger)

  # Get column names.
  cn = strsplit(sub('^#', '', header$header[length(header$header)]), ',')[[1]]

  vars = scanTabix(sanger, param = param)[[1]]
  vars = strsplit(vars, split = '\t')
  vars = matrix(unlist(vars), nrow = length(vars), ncol = length(vars[[1]]), 
                dimnames = list(NULL, cn), byrow = TRUE)
  vars = data.frame(vars)
  vars$pos = as.numeric(vars$pos)
  colnames(vars) = sub('^X', '', colnames(vars))
  
  return(vars)

} # read_sanger_vars()
