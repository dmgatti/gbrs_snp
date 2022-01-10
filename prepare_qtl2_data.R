################################################################################
# Gather founder and sample genotypes and write out in qtl2 format.
# Also gather markers and create json control file.
# Hard-coding the covariates and phenotypes for now.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-05
################################################################################

##### LIBRARIES #####

library(data.table)
library(mmconvert)
library(qtl2convert)
library(qtl2)

source('check_file_access.R')

##### VARIABLES #####

args = commandArgs(trailing = TRUE)

if(length(args) != 3) {

  stop(paste('ERROR: Three arguments required: in_dir, sanger_file, out_path'))

} # if(length(args) != 3)

in_dir      = args[1]
sanger_file = args[2]
out_path    = args[3]

# Test code
#in_dir       = '/fastscratch/dgatti/genotypes'
#sanger_file  = '/projects/compsci/USERS/dgatti/data/gbrs_snp/sanger_transcript_snps_indels_ens102_b38.tsv.bgz'
#out_path     = '/projects/compsci/USERS/dgatti/projects/gbrs_snp/results/qtl2'

check_file_access(in_dir,      4)
check_file_access(sanger_file, 4)
check_file_access(out_path,    2)

# Founder names in Sanger transcript file.
founder_names = c('A_J', 'C57BL_6J', '129S1_SvImJ', 'NOD_ShiLtJ', 
                  'NZO_HlLtJ', 'CAST_EiJ', 'PWK_PhJ', 'WSB_EiJ')

# Covariate file for Attie Islet data.
covar_file = '/projects/compsci/USERS/dgatti/data/gbrs_snp/15NGS-003-gac_DOisletRNAseq_SampleInfo KEL2107-40-35664.csv'

##### PROGRAM #####

# Read in the sample genotypes and merge them together.
files = dir(in_dir, pattern = '_genotypes\\.csv$', full.names = TRUE)

geno = NULL
for(f in files) {

  tmp = fread(f)
  sample = sub('_genotypes\\.csv$', '', basename(f))

  print(sample)
  
  tmp = tmp[,c('marker', 'gt')]
  colnames(tmp)[2] = sample
  
  if(is.null(geno)) {
    geno = tmp
  } else {
    # Union of SNPs.
#    geno = merge(geno, tmp, by = 'marker', all = TRUE, sort = FALSE)
    # Intersection of SNPs.
    geno = merge(geno, tmp, by = 'marker', all = FALSE, sort = FALSE)
  } # else

} # for(f)

# Set genotypes containing '-' to NA.
# TBD: Handle missing genotypes as informative?
for(i in 2:ncol(geno)) {
  geno[[i]][grep('-', geno[[i]])] = NA
}

print(paste(nrow(geno), 'markers in union of genotypes.'))

# Get the founder variants at the same markers.
vars = fread(paste('gunzip -cq', sanger_file), sep = '\t')
# TBD: Use '\t' as column name delimiter in Sanger Tabix file.
cn = strsplit(sub('^#', '', colnames(vars)[1]), ',')[[1]]
colnames(vars) = cn

vars = subset(vars, vars$marker %in% geno$marker)
geno = geno[match(vars$marker, geno$marker),]
stopifnot(all(geno$marker == vars$marker))
vars = as.data.frame(vars)

# Get the marker positions.
pmap = vars[,c('marker', 'chr', 'pos')]

# Get the ref & alternate alleles.
alleles = as.matrix(vars[,c('ref', 'alt')])
rownames(alleles) = vars$marker

# Get the founder genotypes and remove the '/' in the genotypes.
founder_geno = as.matrix(vars[,founder_names])
rownames(founder_geno) = rownames(alleles)
founder_geno = sub('/', '', founder_geno)
rm(vars)

# Encode the genotypes together so that the alleles don't get swapped.
stopifnot(all(geno$marker == rownames(founder_geno)))
tmp = cbind(founder_geno, as.matrix(geno[,-1]))
tmp = encode_geno(geno = tmp, allele_codes = alleles)

founder_geno = data.frame(marker = rownames(tmp), tmp[,1:length(founder_names)])
geno         = data.frame(marker = rownames(tmp), tmp[,(length(founder_names)+1):ncol(tmp)])
rm(tmp)

print('Writing founder & sample genotypes.')

# Write out founder and sample genotypes.
out_file = file.path(out_path, 'founder_geno.csv')
write.csv(founder_geno, file = out_file, quote = FALSE, row.names = FALSE)

out_file = file.path(out_path, 'sample_geno.csv')
write.csv(geno, file = out_file, quote = FALSE, row.names = FALSE)

rm(founder_geno)

# Write out maps.
print('Writing pmap & gmap files.')

pmap$pos = pmap$pos * 1e-6
out_file = file.path(out_path, 'pmap.csv')
write.csv(pmap, file = out_file, quote = FALSE, row.names = FALSE)

# NOTE: mmconvert expects columns in order chr, pos, marker.
tmp = mmconvert(pmap[,c('chr', 'pos', 'marker')], input_type = 'Mbp')
xchr      = which(tmp$chr == 'X')
gmap = tmp[,c('marker', 'chr', 'cM_coxV3_ave')]
colnames(gmap)[3] = 'cM'
gmap$cM[xchr] = tmp[xchr, 'cM_coxV3_female']

out_file = file.path(out_path, 'gmap.csv')
write.csv(gmap, file = out_file, quote = FALSE, row.names = FALSE)

# Write out covariates & phenotypes.
print('Writing covariates & phenotypes.')

meta = read.csv(covar_file)
meta = meta[,c('Customer.Sample.Name', 'NYGC.Sample.ID', 'Sex')]
meta$id = paste('Sample', gsub('_', '.', meta$NYGC.Sample.ID), sep = '_')

# Not sure about generation. Using 25.
covar = data.frame(id  = meta$id, 
                   sex = meta$Sex,
                   gen = 25)
covar = subset(covar, id %in% colnames(geno)[-1])
covar = covar[match(colnames(geno)[-1], covar$id),]
stopifnot(covar$id == colnames(geno)[-1])

out_file = file.path(out_path, 'covar.csv')
write.csv(covar, file = out_file, quote = FALSE, row.names = FALSE)

pheno = data.frame(id   = covar$id,
                   junk = rnorm(nrow(covar)))

out_file = file.path(out_path, 'pheno.csv')
write.csv(pheno, file = out_file, quote = FALSE, row.names = FALSE)

# Write out JSON control file.
print('Writing JSON control file.')

json = '{
  "description": "DO data from Attie Islet",
  "crosstype": "do",
  "sep": ",",
  "na.strings": ["-", "NA"],
  "comment.char": "#",
  "geno": "sample_geno.csv",
  "founder_geno": "founder_geno.csv",
  "gmap": "gmap.csv",
  "pmap": "pmap.csv",
  "pheno": "pheno.csv",
  "covar": "covar.csv",
  "alleles": ["A", "B", "C", "D", "E", "F", "G", "H"],
  "x_chr": "X",
  "genotypes": {
    "A": 1,
    "H": 2,
    "B": 3
  },
  "geno_transposed": true,
  "founder_geno_transposed": true,
  "sex": {
    "covar": "sex",
    "F": "female",
    "M": "male"
  },
  "cross_info": {
    "covar": "gen"
  }
}
'

out_file = file.path(out_path, 'control.json')
writeLines(json, con = out_file)

