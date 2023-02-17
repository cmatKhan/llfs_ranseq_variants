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

pedigree = read_csv("data/triplet_visit2_version3.csv")

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

compare_results = Sys.glob("/mnt/scratch/llfs_rna_dna_compare_test/*/*_match_metrics.csv")

compare_results_df = map(compare_results,read_csv) %>%
  do.call('rbind',.)

# write_rds(compare_results_df, "data/complete_visit1_same_same_20230201.rds")

# write_rds(compare_results_df, "data/intermediate_all_same_compare_20230124.rds")

#compare_results_df = readRDS("data/intermediate_all_same_compare_20230124.rds")

x = compare_results_df %>%
  left_join(distinct(pedigree,subject,gpedid,relative,control),by = c('rna_sample' = 'subject')) %>%
  dplyr::rename(rna_ped = gpedid,rna_related = relative, rna_control = control) %>%
  left_join(distinct(pedigree,subject,gpedid,relative,control),by = c('dna_sample' = 'subject')) %>%
  dplyr::rename(dna_ped = gpedid,dna_related = relative, dna_control = control)

summarized_compare_df = compare_results_df %>%
  mutate(group = paste(rna_sample,rna_visit,dna_sample,sep="_")) %>%
  mutate(group = ifelse(rna_sample == dna_sample, paste0(group,'_same'), paste0(group,'_diff'))) %>%
  group_by(group) %>%
  summarize(total_variants = sum(overlap_fltr),
            total_fltr_match = sum(n_match_fltr),
            total_homo_expr_cand_fltr = sum(homo_expr_cand_fltr)) %>%
  mutate(match_ratio = total_fltr_match/total_variants) %>%
  mutate(same_sample = str_detect(group,"_same"), TRUE, FALSE)

summarized_compare_df %>%
  ggplot(aes(same_sample,match_ratio)) +
  geom_boxplot() +
  ggtitle('chr1-22') +
  coord_cartesian(ylim =  c(0,1))

compare_results_df %>%
  group_by(dna_sample,rna_sample) %>%
  mutate(same_sample = rna_sample==dna_sample) %>%
  ggplot(aes(same_sample,filtered_match_ratio)) + geom_boxplot() +
  ggtitle('chr1-22') +
  coord_cartesian(ylim =  c(0,1))

x = x %>%
  group_by(dna_sample,rna_sample) %>%
  mutate(comparison = ifelse(rna_sample==dna_sample,'Identical',ifelse(rna_ped==dna_ped,'Family','Diff_Family')))

x %>%
  ggplot(aes(comparison,filtered_match_ratio)) + geom_boxplot() +
  ggtitle('chr1-22') +
  coord_cartesian(ylim =  c(0,1))

compare_results_df %>%
  filter(chr=='21') %>%
  mutate(same_sample = rna_sample==dna_sample) %>%
  ggplot(aes(same_sample,filtered_match_ratio)) + geom_boxplot() +
  ggtitle('chr 21 only') +
  coord_cartesian(ylim =  c(0,1))

remaining_comparisons = comparisons_all %>%
  filter(!rna_visit %in% c('11221_1','11221_2',
                          '20686_1','20686_2',
                          '20825_2',
                          '20903_1','20903_2',
                          '21389_2',
                          '25785_1','25785_2',
                          '26715_1',
                          '3751_2',
                          '4263_2',
                          '4307_2')) %>%
  unite(tmp, c('rna_visit','dna_subject'),remove=FALSE)

already_compared = summarized_compare_df %>%
  separate(group, c('rna','visit','dna','label')) %>%
  unite(tmp,c('rna','visit','dna'))

remaining_comparisons %>%
  filter(!tmp %in% already_compared$tmp) %>%
  select(-tmp) %>%
  separate(rna_visit,c('rna_subject','rna_visit')) %>%
  write_tsv("/mnt/scratch/llfs_rna_dna_compare_test/lookups/sex_mislabel_all_comparison_remaining_20220109.tsv",
            col_names = FALSE)
