library(SeqArray)
library(SNPRelate)
library(SeqVarTools)
library(tidyverse)
library(parallel)
# see http://bioconductor.org/packages/release/bioc/vignettes/SeqArray/inst/doc/OverviewSlides.html#(9)
# create a computing cluster with 4 cores
# cl = makeForkCluster(5)
# seqParallelSetup(cl)
# Note: with the pruned set, had error in sample mismatch?
# seqVCF2GDS(
#   "data/pruned_filtered_full_merge.vcf.gz",
#   "data/pruned_filtered_full_merge.gds",
#   reference = 'GRCh38',
#   parallel=TRUE,
#   verbose = TRUE)

# # get variants which have certain properties
# seqGetData(file, '')
#
# # seqSetFilter by variant
#
# # seqSetFilter visit 1
# # extract gentoypes from visit 1
# # seqSetfilter visit 2
# # extract genotypes from visit 2
# # compare and select only variants where visit 1 and visit 2 agree?
#
# # reset filter
# # set filter variants and visit 1 + visit 2 for those mislabeled samples
# #
#
# stopCluster(cl)

file = seqOpen("data/seq_array_attemp2_full_merge.gds")

seqSetFilterChrom(file,include=seq(1,22))

visit1_samples_indicies = which(str_detect(seqGetData(file,'sample.id'), 'visit_1|visit1'))
seqSetFilter(file,sample.sel = visit1_samples_indicies)

subject_df = tibble(tmp = seqGetData(file,'sample.id')) %>%
  separate(tmp,c('id','visit'),sep="_",extra='merge')

variant_depth = seqGetData(file,'annotation/info/DP')
summary(variant_depth)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#2.0      2.0      7.0    463.4     26.0 762562.0
sel_variant_depth = which(variant_depth > 10 & variant_depth < 1000)
summary(variant_depth[sel_variant_depth])

# variant_dist_bias = seqGetData(file, 'annotation/info/VDB')
# summary(variant_dist_bias)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# # 0.00000 0.03496 0.26000 0.34061 0.59801 1.00000
# sel_variant_dist_bias = which(variant_dist_bias > .1)
# summary(variant_dist_bias[sel_variant_dist_bias])

summary(qual(file))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 3.010   7.308  10.792  21.407  26.424 228.432
sel_by_qual = which(qual(file) > 20)

# select only bi-allelic sites
sel_biallelic = which(seqNumAllele(file) == 2)

# select only snv sites
snv_variants = which(isSNV(file))

#sel_variants = intersect(intersect(intersect(sel_variant_depth, sel_variant_dist_bias), sel_by_qual), sel_biallelic)
#sel_variants = intersect(intersect(sel_variant_depth, sel_by_qual), sel_biallelic)
sel_variants = intersect(intersect(intersect(sel_variant_depth,sel_biallelic),snv_variants), sel_by_qual)
seqSetFilter(file, variant.sel = sel_variants)

cl = makeForkCluster(10)
seqParallelSetup(cl)
mac = seqAlleleCount(file,minor=TRUE,parallel = TRUE)
maf = seqAlleleFreq(file,minor=TRUE, parallel = TRUE)
stopCluster(cl)
summary(mac)
summary(maf)

# snpset <-
  # snpgdsLDpruning(
  #   file,
  #   sample.id = seqGetData(file,'sample.id'),
  #   snp.id = seqGetData(file,'variant.id'),
  #   method="corr",
  #   slide.max.bp=10e6,
  #   ld.threshold=sqrt(0.1),
  #   verbose=TRUE,
  #   num.thread=1)

#seqSetFilterCond(file,maf=.3,parallel=TRUE)

het_loc = list()
i = 1
for (r in rownames(dosages)){
  for (c in colnames(dosages)){
    if (dosages[r,c] == 1){
      het_loc[[i]] = c(r,c)
      i = i+1
    }
  }
}

# snpset <-
#   snpgdsLDpruning(
#     file,
#     sample.id = seqGetData(file,'sample.id'),
#     snp.id = seqGetData(file,'variant.id'),
#     method="corr",
#     slide.max.bp=10e6,
#     ld.threshold=sqrt(0.1),
#     verbose=TRUE,
#     num.thread=1)
# this took about 5 hours

#saveRDS(snpset, "data/ld_passing_snp_after_filter.rds")
#snpset = readRDS("data/ld_passing_snp_after_filter.rds")
seqSetFilter(file,variant.sel = unlist(snpset))

x = df %>% group_by(gpedid) %>% group_split() %>% .[[1]]

# try IBD on single family
seqSetFilter(file,sample.sel = c(3,11,148,234))
ibd_king = snpgdsIBDKING(file,
                         sample.id = c("10189066_visit_1",
                                       "10326999_visit_1",
                                       "12075785_visit_1",
                                       "13548632_visit_1"),
                         snp.id = unlist(snpset),
                         num.thread=10)

ibd_mom = snpgdsIBDMoM(file,
                         sample.id = c("10189066_visit_1",
                                       "10326999_visit_1",
                                       "12075785_visit_1",
                                       "13548632_visit_1"),
                         snp.id = unlist(snpset),
                       kinship=TRUE,
                         num.thread=10,
                       verbose = TRUE)
