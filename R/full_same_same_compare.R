library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(SeqVarTools)
library(tidyverse)

llfs_gds = seqOpen('data/LLFS.WGS.freeze5.chr21.gds')
wgs_sample_ids = seqGetData(llfs_gds,'sample.id')
seqClose(llfs_gds)

phantom_rnaseqed = read_csv('data/phantom-rnaseqed.csv')
phantom_donor = read_csv('data/phantoms_donors.csv')

pedigree <- read_csv("data/triplet_visit2_version3.csv")

complete_visit1_same_same <-
  readRDS("data/complete_visit1_same_same_20230201.rds")

sex_mislabel_attempt1 <- readRDS("data/sex_mislabel_all_v_all_attempt1.rds")
sex_mislabel_attempt2 <- readRDS("data/sex_mislabel_all_v_all_attempt2.rds")

s3_vcf_all <- read_csv("data/visit_1_visit_2_filtered_vcf_lookup.txt",
                       col_names = "vcf") %>%
  filter(str_detect(vcf, "tbi$"), negate = TRUE) %>%
  mutate(subject = str_remove(basename(vcf), "_visit\\d_T1.haplotypecaller.filtered.vcf.gz.*"),
         visit = str_extract(vcf, "visit_1|visit_2"),
         in_wgs_dna = ifelse(subject %in% wgs_sample_ids, TRUE, FALSE))

s3_vcf_all_not_in_wgs = filter(s3_vcf_all, !in_wgs_dna)

summarize_comparisons <- function(df) {
  df %>%
    mutate(group = paste(rna_sample, rna_visit, dna_sample, sep = "_")) %>%
    mutate(group = ifelse(rna_sample == dna_sample,
                          paste0(group, "_same"),
                          paste0(group, "_diff"))) %>%
    group_by(group) %>%
    summarize(
      total_variants = sum(overlap_fltr),
      total_fltr_match = sum(n_match_fltr),
      total_homo_expr_cand_fltr = sum(homo_expr_cand_fltr)
    ) %>%
    mutate(match_ratio = total_fltr_match / total_variants) %>%
    mutate(same_sample = str_detect(group, "_same"), TRUE, FALSE)
}
