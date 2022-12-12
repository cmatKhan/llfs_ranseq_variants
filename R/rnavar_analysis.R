library(SeqArray)
library(SNPRelate)
library(SeqVarTools)
library(tidyverse)

# see http://bioconductor.org/packages/release/bioc/vignettes/SeqArray/inst/doc/OverviewSlides.html#(9)
# create a computing cluster with 4 cores
# cl = makeForkCluster(5)
# seqParallelSetup(cl)
# Note: with the pruned set, had error in sample mismatch?
# seqVCF2GDS(
#   "/mnt/scratch/llfs_variant_calling/results/variant_calling/1086_T1/1086_T1.haplotypecaller.filtered.vcf.gz",
#   "data/1086.fltr.gds",
#   reference = 'GRCh38',
#   parallel=TRUE,
#   verbose = TRUE)

file = seqOpen("data/1086.fltr.gds")
seqSetFilterChrom(file,21L)

llfs_gds = seqOpen('data/LLFS.WGS.freeze5.chr21.gds')
seqSetFilter(llfs_gds,sample.id='1404')

# from
getEqualOverlaps = function(vcf1, vcf2) {
  #  This function only makes sense for VCF files with one sample!
  stopifnot(ncol(vcf1) == 1 && ncol(vcf2) == 1)

  #  First take the ranges which are present in both VCF objects
  temp_vcf1 = subsetByOverlaps(vcf1, vcf2)

  #  The indices of each range in vcf2 along the rows of vcf1
  vcf2_indices = match(ranges(temp_vcf1), ranges(vcf2))

  #  Subset those equal ranges to those whose genotype calls match
  temp_vcf1 = temp_vcf1[geno(temp_vcf1)$GT[,1] == geno(vcf2)$GT[vcf2_indices,1],]

  return(temp_vcf1)
}

compare_genotype = function(row){
  seqSetFilterPos(llfs_gds,as.numeric(row[['seqnames']]),as.numeric(row[['start']]))
  seqSetFilterPos(file,as.numeric(row[['seqnames']]),as.numeric(row[['start']]))
  str_replace(as.character(getGenotype(llfs_gds)[1,1]),"\\|","/") == as.character(getGenotype(file)[1,1])
}

file_granges = granges(file)
llfs_granges = granges(llfs_gds)
var_overlap = IRanges::subsetByOverlaps(file_granges,llfs_granges,type='equal')

seqSetFilterPos(file,chr=21L,pos=var_overlap@ranges@start)
seqSetFilterPos(llfs_gds,chr=21L,pos=var_overlap@ranges@start)

get_geno_matrix = function(gds_obj){
  getGenotype(gds_obj)
}

tidy_geno_mat = function(gds_obj,geno_colname){
  getGenotype(gds_obj) %>%
  as_tibble() %>%
  pivot_longer(everything(),values_to=geno_colname,names_to = 'var_id') %>%
  mutate(!!rlang::sym(geno_colname) := str_replace(!!rlang::sym(geno_colname),"\\|","/")) %>%
  left_join(
    granges(gds_obj) %>%
      as_tibble() %>%
      mutate(var_id = names(granges(gds_obj))) %>%
      dplyr::select(var_id,start) %>%
      dplyr::rename(pos=start)) %>%
  dplyr::select(pos,all_of(geno_colname))
}

compare_df = tidy_geno_mat(file,'rnaseq') %>%
  left_join(tidy_geno_mat(llfs_gds,'dnaseq'),by='pos') %>%
  mutate(comp = rnaseq==dnaseq) %>%
  filter(!is.na(comp)) %>%
  group_by(comp) %>%
  tally()


# sum(compare_vector[!is.na(compare_vector)]) / length(compare_vector[!is.na(compare_vector)])
# [1] 0.8335855
