################################################################################
# Create QUILT reference files of SNPs in transcripts.
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
# 2021-12-04
################################################################################

##### LIBRARIES #####

library(AnnotationHub)
library(ensembldb)
library(Rsamtools)
library(VariantAnnotation)

source('check_file_access.R')

##### VARIABLES #####

args = commandArgs(trailing = TRUE)

if(length(args) != 3) {

  stop(paste('ERROR: Three arguments required: sanger_file, ensembl_version, out_path'))

} # if(length(args) != 3)

sanger_file     = args[1]
ensembl_version = as.numeric(args[2])
out_path        = args[3]

# Test code
#sanger_file     = '/projects/omics_share/mouse/GRCm38/genome/variants/snps_indels/rel_2004_v7/mgp_REL2005_snps_indels.vcf.gz'
#ensembl_version = 102
#out_path        = '/projects/compsci/USERS/dgatti/data/gbrs_snp'

check_file_access(sanger_file, 4)
check_file_access(out_path,    2)

if(is.na(ensembl_version)) {

  stop(paste('ERROR: Ensembl version must be an integer. Found:', args[2]))

} # if(is.na(ensembl_version))

# Output filename for SNPs in transcripts.
out_file = file.path(out_path, 'sanger_transcript_snps_indels_ens102_b38.tsv')

##### PROGRAM #####

# Download the Ensembl DB for the requested ensembl version
hub = AnnotationHub()
ensembl = hub[[names(hub)[hub$title == 'Ensembl 102 EnsDb for Mus musculus']]]

# Get transcript locations from ensembl db.
transcripts = transcripts(ensembl)
transcripts = keepStandardChromosomes(transcripts, pruning.mode = 'coarse')

# Get chromosomes.
trans_chr_names = seqlevels(transcripts)

# Get header information from Sanger file.
header = scanVcfHeader(sanger_file)

# Get the chromosome names from the Sanger file.
vcf_chr_names = seqlevels(header)
stopifnot(vcf_chr_names %in% trans_chr_names)

# Sort the chromosomes.
vcf_chr_names = vcf_chr_names[order(as.numeric(vcf_chr_names))]

strains = samples(header)[c(4, 2, 34, 37, 19, 40, 51)]

# Process one chromosome at a time.
for(chr in vcf_chr_names) {

  print(paste('Chromosome', chr))

  # Get all variants on this chromosome.
  gr = GRanges(seqnames = chr, IRanges(start = 0, end = 200e6))
  param = ScanVcfParam(fixed = c('REF', 'ALT'), info = 'INDEL', 
                       samples = strains, geno = 'GT', 
                       which = gr)
  vcf = readVcf(file = sanger_file, param = param)

  print(paste('Found', length(vcf), 'variants.'))

  # For now, keep only SNPs.
  vcf = subset(vcf, info(vcf)$INDEL == FALSE)

  print(paste('Found', length(vcf), 'SNPs'))

  # Intersect variants with transcripts.
  trans_chr = subset(transcripts, seqnames(transcripts) == chr)
  vcf       = subsetByOverlaps(vcf, trans_chr)
  trans_chr = subsetByOverlaps(trans_chr, vcf)

  print(paste(length(vcf), 'variants intersect with', length(trans_chr), 
        'transcripts.'))
        
  # This expands the number of rows where there is more than one alternate
  # alleles. The ExpandedVCF will have one ALT allele per row.
  vcf = expand(vcf, row.names = TRUE)
  
  # Filter to retain polymorphic variants.
  gt  = geno(vcf)$GT
  vcf = subset(vcf, rowMeans(gt == '0/0') < 1.0)
  rm(gt)

  print(paste('Found', length(vcf), 'polymorphic variants.'))

  # QUILT requires one alternate allele. 
  # TBD: Change this to handle mulitple alternate alleles in our own HMM?
  gt      = geno(vcf)$GT
  alleles = apply(gt, 1, unique, simplify = FALSE)
  # NOTE: We only have 7 strains. C57BL/6J is not in the genotypes
  #       because it's the reference. Also, there are 40+ strains in
  #       the file, so we could have quadramorphic SNPs with 1, 2, or 3 
  #       alternate alleles. Other strains may contribute the other 
  #       alleles, but the DO founders might have only 2 alleles. 
  #       So, we could see all '1/1', all '2/2', all 3/3', 
  #       or any of those mixed with '0/0'. I'm going to 
  #       exclude heterozygous calls for now.
  vcf     = subset(vcf, sapply(alleles, function(z) { all(z %in% c('0/0', '1/1', '2/2', '3/3')) }))

  # At this point, we expect two or more homozygous alleles. If all strains are '1/1',
  # then we have one allele in the file, but two in the DO founders because
  # C57BL/6J is '0/0'. The same applies to '2/2' or '3/3'. However, if
  # we have more than two alleles, then we have a non-biallelic SNP and will
  # exclude it (for now). 
  gt      = geno(vcf)$GT
  alleles = apply(gt, 1, unique)
  alleles = lapply(alleles, function(z) { z[z != '0/0'] })
  vcf     = subset(vcf, sapply(alleles, length) == 1)
  # When we expanded the VCF, some rows were duplicated. Remove them.
  gt      = geno(vcf)$GT
  vcf     = subset(vcf, !duplicated(rownames(gt)))
  rm(gt, alleles)

  print(paste('Found', length(vcf), 'biallelic variants.'))
  
  # Verify that we have no more than two homozygous alleles at each SNP.
  gt      = geno(vcf)$GT
  alleles = apply(gt, 1, unique)
  stopifnot(all(unique(unlist(alleles)) %in% c('0/0', '1/1', '2/2', '3/3')))
  stopifnot(all(sapply(alleles, length) <= 2))

  # We need to change the '2/2' and '3/3' to '1/1' since we have an ExpandedVCF
  # with one allele per row. I'm not entirely clear why expand() didn't handle
  # this. It may be because you can still have two alleles in the genotypes.
  gt = gsub(pattern = '^2/2$', replacement = '1/1', gt)
  gt = gsub(pattern = '^3/3$', replacement = '1/1', gt)
  alleles = apply(gt, 1, unique)
  stopifnot(all(unique(unlist(alleles)) %in% c('0/0', '1/1')))
  # Add in C57BL/6J.
  # I can't figure out how to add a column to a VCF object. So I'm recreating it.
  gt = cbind(gt[,1, drop = FALSE], C57BL_6J = '0/0', gt[,-1])
  cd = colData(vcf)
  cd = DataFrame(Samples = 1:8, row.names = c(rownames(cd)[1], 'C57BL_6J', rownames(cd)[-1]))
  
  vcf = VCF(rowRanges = rowRanges(vcf), 
            colData   = cd, 
            exptData  = metadata(vcf),
            fixed     = fixed(vcf),
            info      = info(vcf),
            geno      = SimpleList(GT = gt),
            collapsed = FALSE,
            verbose   = FALSE)

  # Final check for biallelic SNP.
  stopifnot(max(nchar(fixed(vcf)$ALT)) == 1)

  # Write out VCF. We'll need this to make IMPUTE reference files for QUILT.
  writeVcf(vcf, file = file.path(out_path, paste0('sanger_chr', chr, '.vcf')), index = TRUE)
  # Not sure what warning here means yet...

  # Extract positions & genotypes.
  gt  = geno(vcf)$GT
  ref = as.character(fixed(vcf)$REF)
  alt = as.character(fixed(vcf)$ALT)
  for(j in 1:nrow(gt)) {
    gt[j,] = gsub(pattern = '0', replacement = ref[j], gt[j,])
    # We changes all ALT alleles to '1'.
    gt[j,] = gsub(pattern = '1', replacement = alt[j], gt[j,])
  } # for(j)
  gt  = data.frame(gt[,1,drop = FALSE], C57BL_6J = paste(ref, ref, sep = '/'),
                   gt[,2:ncol(gt)])
  for(i in 1:ncol(gt)) {
    gt[[i]] = gsub('\\.', '-', gt[[i]])
  } # for(i)
  indel = info(vcf)$INDEL
  positions = rowRanges(vcf)

  rm(vcf)
  gc()

  # Align transcripts with variants.
  ol = findOverlaps(positions, trans_chr)

  # Transcripts
  tr = trans_chr$tx_id[subjectHits(ol)]
  tr = split(tr, queryHits(ol))
  tr = lapply(tr, unique)
  tr = sapply(tr, paste0, collapse = ';')

  # Genes
  ge = trans_chr$gene_id[subjectHits(ol)]
  ge = split(ge, queryHits(ol))
  ge = lapply(ge, unique)
  ge = sapply(ge, paste0, collapse = ';')

  # Transcript strand
  st = CharacterList(strand(trans_chr))[[1]]
  st = st[subjectHits(ol)]
  st = split(st, queryHits(ol))
  st = lapply(st, unique)
  st = sapply(st, paste0, collapse = ';')

  # Assemble the output.
  output = data.frame(marker = paste0(chr, ':', start(positions), '_', ref, '/', alt),
                      chr    = seqnames(positions),
                      pos    = start(positions),
                      ref    = ref,
                      alt    = alt,
                      tx_id  = tr,
                      tx_strand = st,
                      gene_id = ge,
                      indel  = indel,
                      gt
                     )
  # Fix the 129S1 column name.
  colnames(output) = sub('^X', '', colnames(output))

  output = output[order(output$pos),]

  print(paste('Writing', nrow(output), 'variants.'))

  # Write out a header for the file if this is chr 1.
  if(chr == '1') {

    h = rep('', 8)
    h[1] = '# Sanger SNPs & Indels in Transcripts'
    h[2] = '# Genome: GRCm38'
    h[3] = '# Sanger: version 7'
    h[4] = '# Ensembl: 102'
    h[5] = '# Creator:Daniel Gatti'
    h[6] = '# Email: dan.gatti@jax.org'
    h[7] = '# Date: 2021-12-04'
    h[8] = paste0('#', paste0(colnames(output), collapse = ','))

    writeLines(text = h, con = out_file)

  } # if(chr == '1')

  write.table(output, file = out_file, sep = '\t', quote = FALSE, row.names = FALSE, 
              col.names = FALSE, append = TRUE)

  rm(positions, gt, ref, alt, tr, ge, st, indel, output)
  gc()

} # for(chr)

# Zip and index the tabix file.
print(paste('Zipping and indexing Tabix file.'))

zipfile = bgzip(out_file, overwrite = TRUE)
indexTabix(file = zipfile, seq = 2L, start = 3L, end = 3L, skip = 8L)

