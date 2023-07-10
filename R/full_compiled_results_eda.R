library(tidyverse)

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

# Now done via script on the cluster
#
# gather_data_from_scratch = function(){
#   compare_results = Sys.glob("/mnt/scratch/llfs_rna_dna_compare_test/tmp_vcf_gds/*/*_match_metrics.csv")
#
#   map(compare_results,read_csv) %>%
#     do.call('rbind',.)
# }
#
# gather_data_from_genoarr = function(){
#   compare_results = Sys.glob("/mnt/scratch/llfs_rna_dna_compare_test/genoarr_comp/*/*_match_metrics.csv")
#
#   map(compare_results,read_csv) %>%
#     do.call('rbind',.)
# }


complete_visit1_same_same <-
  readRDS("data/complete_visit1_same_same_20230201.rds")

sex_mislabel_attempt1 <-
  readRDS("data/sex_mislabel_all_v_all_attempt1.rds")
sex_mislabel_attempt2 <-
  readRDS("data/sex_mislabel_all_v_all_attempt2.rds")

pedigree = read_csv("data/triplet_visit2_version3.csv") %>%
  mutate(subject = as.character(subject))

sex_mislabels = read_csv('data/suspicious_sex_samples.csv')

#full_same_same_comp = gather_data_from_scratch()
#write_csv(full_same_same_comp,"data/full_same_same_comp_20230218.csv")
full_same_same_comp = read_csv("data/full_same_same_comp_20230220.csv")

# note -- this is here b/c I needed to do an all by all on the clearly
# mislabelled samples from the new batches
full_same_same_comp = read_csv('data/same_same_new_batches_20230313.csv')


# NOTE: must run the full_same_same_compare script for this. it is here b/c
# of badly needed and only partially begun refactoring
parsed_s3 = list(
  lookup = s3_vcf_all
)

all_by_all_newsamples = read_csv('data/all_by_all_newsamples.csv')

summary_all_by_all_newsamples = summarize_comparisons(all_by_all_newsamples) %>%
  separate(group, c('rna_subject', 'visit_tmp', 'visit', 'dna_subject', 'tmp'), sep="_") %>%
  select(-c(same_sample, `TRUE`, `FALSE`, visit_tmp, tmp))

summary_all_by_all_newsamples %>%
  filter(match_ratio > .9) %>%
  group_by(rna_subject) %>%
  filter(n()==1) %>%
  View()

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

match_summary %>%
  group_by(match_level) %>%
  tally()

match_summary %>%
  ggplot(aes(match_level, match_ratio)) +
  geom_boxplot()

# this was used originally for the full set
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

# this is used for the new samples afte rlooking at the boxplots -- there is a
# set below .5 that obviously needs to be checked
all_by_all_candidates =
  filter(match_summary, match_level == 'mislabel') %>%
  mutate(sex_mislabel = subject %in% unique(sex_mislabels$subject)) %>%
  filter(!sex_mislabel) %>%
  select(subject,visit) %>%
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

create_all_by_all_lookup = function(row, wgs_samples_list){

  tibble(rna_subject = row[['subject']],
         visit=row[['visit']],
         rna_gds = row[['rna_gds']],
         wgs_dna = wgs_samples_list)

}

# write_tsv(tibble(wgs = wgs_sample_ids),
#           "/mnt/scratch/llfs_rna_dna_compare_test/all_by_all_newsamples/lookups/dna_samples.txt",
#           col_names = FALSE)
#
# write_tsv(all_by_all_candidates,
#           "/mnt/scratch/llfs_rna_dna_compare_test/all_by_all_newsamples/lookups/rna_samples.txt",
#           col_names = FALSE)

# all_by_all_subj_lookup = apply(all_by_all_candidates,1,create_all_by_all_lookup, wgs_in_rna) %>%
#   do.call('rbind',.) %>%
#   mutate(visit = str_remove(visit,'visit_'))
#
# all_by_all_subj_lookup_split = all_by_all_subj_lookup %>%
#   group_split(grp = as.integer(gl(n(), 10000, n())), .keep = FALSE)
#
# names(all_by_all_subj_lookup_split) =
#   as.character(seq(1,length(all_by_all_subj_lookup_split)))

# gather_data_from_all_by_all = function(){
#   compare_results = Sys.glob("/mnt/scratch/llfs_rna_dna_compare_test/all_by_all_wgs/*/*_match_metrics.csv")
#
#   map(compare_results,read_csv) %>%
#     do.call('rbind',.)
# }

#all_by_all_metrics = gather_data_from_all_by_all()
#all_by_all_metrics %>% write_csv("data/all_by_all_metrics_attempt1_20230221.csv")
all_by_all_metrics = read_csv("data/all_by_all_metrics_rnaseqed_only.csv")

all_by_all_summary = all_by_all_metrics %>%
  dplyr::rename(rna_sample=rna_subject,
                rna_visit=visit,
                dna_sample=wgs_dna) %>%
  summarize_comparisons()

all_by_all_summary = all_by_all_summary %>%
  separate(group, c('rna_subject', 'visit', 'dna_subject', 'tmp')) %>%
  select(-tmp)

best_match = all_by_all_summary %>%
  group_by(rna_subject, visit) %>%
  filter(rank(match_ratio) == n()) %>%
  filter(match_ratio > .9) %>%
  mutate(subject_visit = paste(rna_subject, visit, sep="_"))

all_by_all_reident = best_match %>%
  filter(!same_sample)

View(all_by_all_reident)

sex_mislabel_comb = read_csv("data/sex_mislabel_all_v_all_combined.rds")

sex_mislabel_reident = sex_mislabel_comb %>%
  summarize_comparisons() %>%
  separate(group, c('rna_subject', 'visit', 'dna_subject', 'tmp')) %>%
  select(-tmp) %>%
  group_by(rna_subject, visit) %>%
  filter(rank(match_ratio) == n()) %>%
  filter(match_ratio > .9)

view(sex_mislabel_reident)

reidentified_samples = rbind(all_by_all_reident, sex_mislabel_reident) %>%
  mutate(subject_visit = paste(rna_subject, visit, sep="_"))

# write_csv(reidentified_samples, "data/reidentified_samples_20230305.csv")

full_confirmed_samples = rbind(best_match, sex_mislabel_reident) %>%
  mutate(subject_visit = paste(rna_subject, visit, sep="_"))

rna_samples <- read_csv("data/20221212_rnaseq_metadata.csv") %>%
  mutate(visit = ifelse(visit == 'visit1', 'visit_1', visit)) %>%
  mutate(visit = ifelse(visit == 'visit2', 'visit_2', visit))

wgs_not_in_rnaseq = wgs_sample_ids[!wgs_sample_ids %in% unique(rna_samples$subject)]

match_against_all_wgs = match_summary %>%
  mutate(subject_visit = paste(subject, str_extract(visit,'\\d'), sep="_")) %>%
  filter(!subject_visit %in% unique(full_confirmed_samples$subject_visit)) %>%
  filter(match_level != 'positive')

all_by_all_not_in_rnaseq_candidates = all_by_all_candidates %>%
  right_join(match_against_all_wgs) %>%
  select(subject,visit,rna_gds) %>%
  left_join(select(parsed_s3$lookup, subject,visit,vcf)) %>%
  mutate(vcf = str_replace(vcf, ".vcf.gz", '.gds')) %>%
  mutate(rna_gds = ifelse(is.na(rna_gds), vcf, rna_gds)) %>%
  select(-vcf)

intermediate_compiled = read_csv("/mnt/scratch/llfs_rna_dna_compare_test/intermediate_compiled_metrics.csv")

intermediate_compiled_summary = intermediate_compiled %>%
  summarize_comparisons()

remove_from_lookups = intermediate_compiled %>%
  group_by(rna_sample,rna_visit) %>%
  tally() %>%
  filter(n == 3121) %>%
  mutate(done = TRUE) %>%
  dplyr::select(rna_sample,rna_visit,done) %>%
  dplyr::rename(subject=rna_sample,visit=rna_visit) %>%
  mutate(subject = as.character(subject), visit = as.character(visit))

all_by_all_not_in_rnaseq_candidates %>%
  left_join(remove_from_lookups) %>%
  replace_na(list(done=FALSE)) %>%
  filter(!done) %>%
  dplyr::select(-done)
# %>%
#   write_tsv("/mnt/scratch/llfs_rna_dna_compare_test/all_by_all_not_rnaseqed/lookups/rnaseq_samples_fltr.txt",
#             col_names = FALSE)

# write_tsv(all_by_all_not_in_rnaseq_candidates,
#           "/mnt/scratch/llfs_rna_dna_compare_test/all_by_all_not_rnaseqed/lookups/rnaseq_samples.txt",
#           col_names = FALSE)

# wgs_not_in_rnaseq_df = tibble(wgs = wgs_not_in_rnaseq) %>%
#   write_tsv("/mnt/scratch/llfs_rna_dna_compare_test/all_by_all_not_rnaseqed/lookups/dna_samples.txt",
#             col_names = FALSE)

intermediate_compiled = read_csv("/mnt/scratch/llfs_rna_dna_compare_test/all_by_all_not_rnaseqed_compiled.csv")

intermediate_compiled_summary = intermediate_compiled %>%
  summarize_comparisons()

