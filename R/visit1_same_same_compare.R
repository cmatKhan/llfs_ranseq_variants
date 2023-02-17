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

suspicious_sex_samples = read_csv('data/suspicious_sex_samples.csv')

sex_mislabel_all_v_all = rbind(
  readRDS("data/sex_mislabel_all_v_all_attempt1.rds"),
  readRDS("data/sex_mislabel_all_v_all_attempt2.rds")) %>%
  mutate(rna_dna = paste(rna_sample,dna_sample, sep="_")) %>%
  filter(chr == 1)

visit1_same_same = readRDS("data/complete_visit1_same_same_20230201.rds")

# NOTE: HERE WE'RE LOOKING FOR MATCHES, SO FILTER IS MATCH_RATIO > .9
summarized_sex_mislabel = sex_mislabel_all_v_all %>%
  mutate(group = paste(rna_sample,rna_visit,dna_sample,sep="_")) %>%
  mutate(group = ifelse(rna_sample == dna_sample, paste0(group,'_same'), paste0(group,'_diff'))) %>%
  group_by(group) %>%
  summarize(total_variants = sum(overlap_fltr),
            total_fltr_match = sum(n_match_fltr),
            total_homo_expr_cand_fltr = sum(homo_expr_cand_fltr)) %>%
  mutate(match_ratio = total_fltr_match/total_variants) %>%
  mutate(same_sample = str_detect(group,"_same"), TRUE, FALSE) %>%
  filter(match_ratio > .9)

# NOTE: HERE WE'RE LOOKING FOR MISMATCHES, SO MATCH_RATIO IS < .9
summarized_visit1_same_same = visit1_same_same %>%
  mutate(group = paste(rna_sample,rna_visit,dna_sample,sep="_")) %>%
  mutate(group = ifelse(rna_sample == dna_sample, paste0(group,'_same'), paste0(group,'_diff'))) %>%
  group_by(group) %>%
  summarize(total_variants = sum(overlap_fltr),
            total_fltr_match = sum(n_match_fltr),
            total_homo_expr_cand_fltr = sum(homo_expr_cand_fltr)) %>%
  mutate(match_ratio = total_fltr_match/total_variants) %>%
  mutate(same_sample = str_detect(group,"_same"), TRUE, FALSE) %>%
  filter(match_ratio < .9) %>%
  separate(group, c('rna_sample', 'visit', 'dna_sample', 'remove_me'),remove = FALSE)

summarized_visit1_same_same = summarized_visit1_same_same %>%
  mutate(previously_known =
           ifelse(rna_sample %in% suspicious_sex_samples$subject,
                  TRUE, FALSE)) %>%
  filter(!previously_known) %>%
  View()




