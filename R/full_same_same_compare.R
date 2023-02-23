library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(SeqVarTools)
library(tidyverse)
library(RSQLite)
library(readxl)

llfs_gds <- seqOpen("data/LLFS.WGS.freeze5.chr21.gds")
wgs_sample_ids <- seqGetData(llfs_gds, "sample.id")
seqClose(llfs_gds)

generate_same_same_lookup = function(
    write_scratch_lookup = FALSE,
    write_local_missing_dna_data = FALSE){

  llfs_gds <- seqOpen("data/LLFS.WGS.freeze5.chr21.gds")
  wgs_sample_ids <- seqGetData(llfs_gds, "sample.id")
  seqClose(llfs_gds)

  geno_arr_gds <- seqOpen("/mnt/ref/data/llfs/geno_chip/llfs_gwas.chr22.gds")
  geno_arr_sample_id <- seqGetData(geno_arr_gds, "sample.id")
  seqClose(geno_arr_gds)

  whatdatall <- read_csv("data/whatdatall.csv") %>%
    select(id, subject) %>%
    distinct(subject, .keep_all = TRUE) %>%
    mutate(
      id = as.character(id),
      subject = as.character(subject)
    ) %>%
    dplyr::rename(
      whatdatall_id = id,
      whatdatall_subj = subject
    ) %>%
    filter(complete.cases(.))

  phantom_rnaseqed <- read_csv("data/phantom-rnaseqed.csv") %>%
    select(id, phmid, subject) %>%
    distinct(id, phmid, subject) %>%
    mutate(
      id = as.character(id),
      subject = as.character(subject),
      phmid = as.character(phmid)
    ) %>%
    dplyr::rename(phan_subj = subject) %>%
    select(id, phan_subj, phmid) %>%
    mutate(missing_subject = ifelse(is.na(phan_subj), TRUE, FALSE)) %>%
    left_join(whatdatall, by = c("id" = "whatdatall_id"))

  phantom_rnaseqed %>%
    filter(is.na(phan_subj)) %>%
    left_join(whatdatall, by = c("whatdatall_subj" = "whatdatall_subj"))

  phantom_rnaseqed <- phantom_rnaseqed %>%
    mutate(phan_subj = ifelse(is.na(phan_subj), whatdatall_subj, phan_subj))

  legal_fails <- read_excel(
    file.path(
      "/mnt/lts/personal/chasem/llfs",
      "data_processing/data/LLFS_DKsamples_NOT_RNASeq_20220523.xlsx"
    )
  )

  s3_vcf_all <- read_csv("data/visit_1_visit_2_filtered_vcf_lookup.txt",
    col_names = "vcf"
  ) %>%
    filter(str_detect(vcf, "tbi$", negate = TRUE)) %>%
    mutate(
      subject = str_remove(
        basename(vcf),
        "_visit\\d_T1.haplotypecaller.filtered.vcf.gz.*"
      ),
      visit = str_extract(vcf, "visit_1|visit_2"),
      subject = as.character(subject),
      phantom = ifelse(subject %in%
        as.character(na.omit(phantom_rnaseqed$phmid)),
      TRUE, FALSE
      )
    ) %>%
    left_join(phantom_rnaseqed, by = c("subject" = "phmid")) %>%
    mutate(subject = ifelse(phantom &
      !is.na(phan_subj), phan_subj, subject)) %>%
    mutate(
      in_wgs_dna = ifelse((subject %in% wgs_sample_ids |
        whatdatall_subj %in% wgs_sample_ids),
      TRUE, FALSE
      ),
      in_genoarr_dna = ifelse((subject %in% geno_arr_sample_id |
        whatdatall_subj %in% geno_arr_sample_id),
      TRUE, FALSE)) %>%
    mutate(legal_fail =
             ifelse(subject %in%
                      legal_fails$subject |
                      whatdatall_subj %in% legal_fails$subject,
                    TRUE, FALSE))

  s3_vcf_all_not_in_dna <- filter(s3_vcf_all, !in_wgs_dna, !in_genoarr_dna)

  chr <- "21"
  chr1_wgs_gds <- paste0("/ref/mblab/data/llfs/agds/LLFS.WGS.freeze5.chr", chr, ".gds")
  chr1_arr_gds <- paste0("/ref/mblab/data/llfs/geno_chip/llfs_gwas.chr", chr, ".gds")

  full_same_same_comp <- s3_vcf_all %>%
    filter(in_wgs_dna | in_genoarr_dna) %>%
    mutate(
      gds = ifelse(in_wgs_dna, chr1_wgs_gds, chr1_arr_gds),
      vcf_index = paste0(vcf, ".tbi"),
      chr = chr
    ) %>%
    select(vcf, vcf_index, subject, visit, chr, gds)

  if(write_scratch_lookup){
    write_tsv(
      full_same_same_comp,
      "/mnt/scratch/llfs_rna_dna_compare_test/lookups/full_same_same_comp.txt",
      col_names = FALSE)
  }

  if(write_local_missing_dna_data){
    s3_vcf_all_not_in_dna %>%
      select(-c(vcf, id, phan_subj, missing_subject, whatdatall_subj)) %>%
      write_csv("data/rna_subjects_not_in_dna_20230217.csv")
  }

  list(
    s3_vcf_all = s3_vcf_all,
    lookup = full_same_same_comp
  )
}

#' after collecting the data from scratch, summarize it
summarize_comparisons <- function(df) {
  df %>%
    mutate(group = paste(rna_sample, rna_visit, dna_sample, sep = "_")) %>%
    mutate(group = ifelse(rna_sample == dna_sample,
      paste0(group, "_same"),
      paste0(group, "_diff")
    )) %>%
    group_by(group) %>%
    summarize(
      total_variants = sum(overlap_fltr),
      total_fltr_match = sum(n_match_fltr),
      total_homo_expr_cand_fltr = sum(homo_expr_cand_fltr)
    ) %>%
    mutate(match_ratio = total_fltr_match / total_variants) %>%
    mutate(same_sample = str_detect(group, "_same"), TRUE, FALSE)
}

gather_data_from_scratch = function(){
  compare_results = Sys.glob("/mnt/scratch/llfs_rna_dna_compare_test/tmp_vcf_gds/*/*_match_metrics.csv")

  map(compare_results,read_csv) %>%
    do.call('rbind',.)
}

gather_data_from_genoarr = function(){
  compare_results = Sys.glob("/mnt/scratch/llfs_rna_dna_compare_test/genoarr_comp/*/*_match_metrics.csv")

  map(compare_results,read_csv) %>%
    do.call('rbind',.)
}


complete_visit1_same_same <-
  readRDS("data/complete_visit1_same_same_20230201.rds")

sex_mislabel_attempt1 <- readRDS("data/sex_mislabel_all_v_all_attempt1.rds")
sex_mislabel_attempt2 <- readRDS("data/sex_mislabel_all_v_all_attempt2.rds")

pedigree = read_csv("data/triplet_visit2_version3.csv") %>%
  mutate(subject = as.character(subject))

sex_mislabels = read_csv('data/suspicious_sex_samples.csv')

parsed_s3 = generate_same_same_lookup()

#full_same_same_comp = gather_data_from_scratch()
#write_csv(full_same_same_comp,"data/full_same_same_comp_20230218.csv")
full_same_same_comp = read_csv("data/full_same_same_comp_20230220.csv")

full_same_same_genoarr = gather_data_from_genoarr()

summary_full_same_same_genoarr = summarize_comparisons(full_same_same_genoarr)%>%
  select(-c(same_sample, `TRUE`, `FALSE`))

summary_full_same_same_comp = summarize_comparisons(full_same_same_comp) %>%
  select(-c(same_sample, `TRUE`, `FALSE`))

summary(summary_full_same_same_comp$match_ratio)

summary_full_same_same_comp = summary_full_same_same_comp %>%
  replace(is.na(.),0)

summary(summary_full_same_same_comp$match_ratio)

visit1_visit2_mismatch = summary_full_same_same_comp %>%
  mutate(visit = str_extract(group, 'visit_\\d')) %>%
  mutate(subject = str_extract(group, '^\\d+')) %>%
  select(visit,subject,match_ratio) %>%
  pivot_wider(names_from = visit, values_from = match_ratio) %>%
  mutate(visit_diff = abs(visit_1-visit_2)) %>%
  filter(visit_diff > .01) %>%
  arrange(desc(visit_diff))

match_summary = summary_full_same_same_comp  %>%
  mutate(visit = str_extract(group, 'visit_\\d')) %>%
  mutate(subject = str_extract(group, '^\\d+')) %>%
  select(visit,subject,match_ratio) %>%
  mutate(positive = match_ratio > .97,
         likely = between(match_ratio,.9,.97),
         close = (match_ratio>.8 & match_ratio < .9),
         mislabel = match_ratio<=.8) %>%
  pivot_longer(c(positive,likely,close,mislabel),
               names_to = 'match_level', values_to = 'tmp') %>%
  # filter only those that are TRUE in the match_level column
  filter(tmp) %>%
  select(-tmp) %>%
  left_join(pedigree, by = 'subject') %>%
  mutate(match_level = factor(match_level,
                              levels = c('mislabel','close','likely','positive'),
                              ordered = TRUE))

match_summary %>% group_by(match_level,Twinstatus) %>% tally()

match_summary %>%
  ggplot(aes(match_level, match_ratio)) +
  geom_boxplot()

all_by_all_candidates =
  filter(match_summary, match_level != 'positive') %>%
  mutate(sex_mislabel = subject %in% unique(sex_mislabels$subject)) %>%
  filter(!sex_mislabel) %>%
  select(subject,visit) %>%
  # add a mislabel from the genoarr data
  rbind(tibble(subject = rep('14547',2), visit = paste0('visit_',seq(1,2)))) %>%
  left_join(parsed_s3$lookup) %>%
  select(subject,visit,vcf) %>%
  dplyr::rename(rna_gds = vcf) %>%
  mutate(rna_gds = str_replace(rna_gds,'.vcf.gz','.gds'))

rna_metadata = read_csv("data/20221212_rnaseq_metadata.csv") %>%
  mutate(visit = ifelse(visit == 'visit1', 'visit_1', visit)) %>%
  mutate(visit = ifelse(visit == 'visit2', 'visit_2', visit))

wgs_in_rna = rna_metadata %>%
  filter(subject %in% wgs_sample_ids) %>%
  distinct(subject) %>%
  pull(subject)

create_all_by_all_lookup = function(row){

  tibble(rna_subject = row[['subject']],
         visit=row[['visit']],
         rna_gds = row[['rna_gds']],
         wgs_dna = wgs_in_rna)

}

all_by_all_subj_lookup = apply(all_by_all_candidates,1,create_all_by_all_lookup) %>%
  do.call('rbind',.) %>%
  mutate(visit = str_remove(visit,'visit_'))

all_by_all_subj_lookup_split = all_by_all_subj_lookup %>%
  group_split(grp = as.integer(gl(n(), 10000, n())), .keep = FALSE)

names(all_by_all_subj_lookup_split) =
  as.character(seq(1,length(all_by_all_subj_lookup_split)))

gather_data_from_all_by_all = function(){
  compare_results = Sys.glob("/mnt/scratch/llfs_rna_dna_compare_test/all_by_all_wgs/*/*_match_metrics.csv")

  map(compare_results,read_csv) %>%
    do.call('rbind',.)
}

#all_by_all_metrics = gather_data_from_all_by_all()
#all_by_all_metrics %>% write_csv("data/all_by_all_metrics_attempt1_20230221.csv")
all_by_all_metrics = read_csv("data/all_by_all_metrics_attempt1_20230221.csv") %>%
  dplyr::rename(visit = rna_visit,
                rna_subject = rna_sample,
                wgs_dna = dna_sample) %>%
  mutate(visit = as.character(visit),
         rna_subject = as.character(rna_subject),
         wgs_dna = as.character(wgs_dna))

x = all_by_all_subj_lookup %>%
  left_join(all_by_all_metrics)

y = x %>%
  filter(!complete.cases(.)) %>%
  group_by(rna_subject,visit) %>%
  group_split()

name_y = function(df){
  df %>%
    distinct(rna_subject,visit) %>%
    as.character() %>%
    paste(collapse='_') %>%
    paste('rna',.,sep="_")
}

names(y) = unlist(map(y,name_y))

y = map(y,~distinct(.,wgs_dna))

lookups_output = "/mnt/scratch/llfs_rna_dna_compare_test/all_by_all_wgs/lookups"
map(names(y), ~write_tsv(y[[.]],
                         file.path(lookups_output,paste0(.,'.txt')),
                         col_names = FALSE))

z = x %>%
  filter(!complete.cases(.)) %>%
  group_by(rna_subject,visit) %>%
  group_split() %>%
  map(distinct,rna_subject,visit,rna_gds) %>%
  do.call('rbind',.) %>%
  write_tsv(file.path(lookups_output,"rna_subj_lookup.txt"),col_names = FALSE)

# map(names(all_by_all_subj_lookup_split),
#     ~write_tsv(all_by_all_subj_lookup_split[[.]],
#                paste0("/mnt/scratch/llfs_rna_dna_compare_test/",
#                       "all_by_all_wgs/lookups/subject_lookup_",
#                       as.character(.),
#                       ".txt"),
#             col_names = FALSE))
#
# chunk9 = read_csv("/mnt/scratch/llfs_variant_calling/samplesheet/completed/visit_2_chunk_9_samplesheet.csv") %>%
#   separate(sample, c('subject','visit')) %>%
#   mutate(visit = 'visit_2')
#
# x = parsed_s3$lookup %>%
#   right_join(chunk9)
#
# full_same_same_comp %>%
#   dplyr::rename(subject = rna_sample, visit = rna_visit) %>%
#   mutate(subject = as.character(subject),
#          chr = as.character(chr)) %>%
#   full_join(rna_metadata) %>%
#   filter(is.na(chr)) %>%
#   filter(visit != 'control') %>%
#   left_join(parsed_s3$s3_vcf_all) %>%
#   filter(in_genoarr_dna) %>%
#   mutate(vcf_index = paste0(vcf,".tbi")) %>%
#   select(vcf,vcf_index, subject,visit) %>%
#   write_tsv("/mnt/scratch/llfs_rna_dna_compare_test/genoarr_comp/subject_lookup.txt")


















