library(tidyverse)

#' Summarize Comparisons
#'
#' Summarizes various statistics for each group of RNA and DNA samples.
#'
#' @param df A data frame containing columns 'rna_sample', 'rna_visit',
#'           'dna_sample', 'overlap_fltr', 'n_match_fltr', and
#'           'homo_expr_cand_fltr'.
#'
#' @return A data frame containing columns 'group', 'total_variants',
#'         'total_fltr_match', 'total_homo_expr_cand_fltr', 'match_ratio',
#'         and 'same_sample'.
#'
#' @details The 'group' column is created by concatenating the values of
#' 'rna_sample', 'rna_visit', and 'dna_sample' with underscores. If 'rna_sample'
#' and 'dna_sample' are equal, '_same' is appended to the end of the string,
#' otherwise '_diff' is appended. The function then groups the data frame by
#' 'group' and calculates the sum of 'overlap_fltr', 'n_match_fltr', and
#' 'homo_expr_cand_fltr' for each group. The 'match_ratio' column is calculated
#' as the ratio of 'total_fltr_match' to 'total_variants'. Finally, the
#' 'same_sample' column is created using str_detect to check if the 'group'
#' column contains the string '_same', returning TRUE if it does and FALSE
#' otherwise.
#'
#' @examples
#' df <- data.frame(
#'   rna_sample = c("RNA1", "RNA2", "RNA3", "RNA4"),
#'   rna_visit = c("V1", "V2", "V3", "V4"),
#'   dna_sample = c("DNA1", "DNA2", "DNA3", "DNA4"),
#'   overlap_fltr = c(10, 20, 30, 40),
#'   n_match_fltr = c(2, 4, 6, 8),
#'   homo_expr_cand_fltr = c(1, 3, 5, 7)
#' )
#' summarize_comparisons(df)
#'
#' @import dplyr
#' @import stringr
summarize_comparisons <- function(df) {
  df %>%
    mutate(group = paste(rna_sample, rna_visit, dna_sample, sep = "_")) %>%
    mutate(group = ifelse(rna_sample == dna_sample,
      paste0(group, "_same_", row_number()),
      paste0(group, "_diff_", row_number())
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

# load data -------------------------------------------------------------------

pedigree = read_csv("data/triplet_visit2_version3.csv") %>%
  mutate(subject = as.character(subject))

# compile full comparison results ---------------------------------------------

all_by_all_metrics_rnaseqed <-
  read_csv("data/all_by_all_metrics_rnaseqed_only.csv")

sex_mislabel_comb <- read_csv("data/sex_mislabel_all_v_all_combined.rds")
colnames(sex_mislabel_comb) <- colnames(all_by_all_metrics_rnaseqed)

intermediate_path <- file.path(
  "/mnt/scratch/llfs_rna_dna_compare_test",
  "intermediate_compiled_metrics.csv"
)
intermediate_compiled <- read_csv(intermediate_path)
colnames(intermediate_compiled) <- colnames(all_by_all_metrics_rnaseqed)

all_by_all_path_not_rnaseqed_path <- file.path(
  "/mnt/scratch/llfs_rna_dna_compare_test",
  "all_by_all_not_rnaseqed_compiled.csv"
)
all_by_all_not_rnaseqed <- read_csv(all_by_all_path_not_rnaseqed_path)
colnames(all_by_all_not_rnaseqed) <- colnames(all_by_all_metrics_rnaseqed)

full_res <- rbind(
  all_by_all_metrics_rnaseqed,
  sex_mislabel_comb,
  intermediate_compiled,
  all_by_all_not_rnaseqed) %>%
  unite('subject_visit',c(rna_subject, visit),  remove=FALSE) %>%
  group_by(rna_subject, visit, wgs_dna) %>%
  distinct(.keep_all = TRUE) %>%
  dplyr::rename(rna_sample=rna_subject,
                rna_visit=visit,
                dna_sample=wgs_dna)

summary_full_res = summarize_comparisons(full_res)

relabelled_samples_tranch1 = summary_full_res %>%
  mutate(group = str_remove(group, "visit_")) %>%
  separate_wider_delim(group, names=c('subject', 'visit',
                                      'dna_relabel', 'tmp'),
                       delim="_", too_many='merge') %>%
  group_by(subject, visit) %>%
  summarize(max_match_ratio = max(match_ratio),
            num_match_ratio_greater_90 = sum(match_ratio>.9),
            samples_greater_90 = paste(dna_relabel[match_ratio > 0.9], collapse = ",")) %>%
  ungroup()

all_by_all_newsamples = 'data/all_by_all_newsamples.csv'
all_by_all_newsamples <- read_csv(all_by_all_newsamples)

summary_all_by_all_newsamples = summarize_comparisons(all_by_all_newsamples)

relabelled_samples_tranch2 = summary_all_by_all_newsamples %>%
  mutate(group = str_remove(group, "visit_")) %>%
  separate_wider_delim(group, names=c('subject', 'visit',
                                      'dna_relabel', 'tmp'),
                       delim="_", too_many='merge') %>%
  group_by(subject, visit) %>%
  summarize(max_match_ratio = max(match_ratio),
            num_match_ratio_greater_90 = sum(match_ratio>.9),
            samples_greater_90 = paste(dna_relabel[match_ratio > 0.9], collapse = ",")) %>%
  ungroup()
###############################################################################

# compile same-same comparison results ----------------------------------------
same_same_list = list(
  tranch1 =
    read_csv("data/full_same_same_comp_20230220.csv"),
  tranch2 =
    read_csv("data/same_same_new_batches_20230313.csv")
)

# examine the same same results -----------------------------------------------

same_same_summary = map(same_same_list, summarize_comparisons)

# examine the unduplicated set
label_matches = function(df){

  df %>%
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
                                ordered = TRUE)) %>%
    mutate(visit = str_extract(visit,'\\d'))

}

same_same_labelled = map(same_same_summary, label_matches)

same_same_labelled$tranch1 %>%
  group_by(match_level) %>%
  tally()

same_same_labelled$tranch2 %>%
  group_by(match_level) %>%
  tally()

same_same_labelled$tranch1 %>%
  ggplot(aes(match_level, match_ratio)) +
  geom_boxplot() +
  ggtitle('tranch1')

same_same_labelled$tranch2 %>%
  ggplot(aes(match_level, match_ratio)) +
  geom_boxplot() +
  ggtitle('tranch2')

# note: all likely matches match themselves within a hundredth or so
# of another sample -- these should be considered correctly labelled.
# only mislabel and close are considered here, therefore
tranch1_relabel_map = same_same_labelled$tranch1 %>%
  filter(match_level %in% c('mislabel', 'close')) %>%
  select(visit,subject,match_ratio) %>%
  dplyr::rename(same_same = match_ratio) %>%
  left_join(relabelled_samples_tranch1,
            by = c('subject', 'visit'),
            multiple='all') %>%
  mutate(dna_relabel = samples_greater_90,
         notes = "")

# tranch1_relabel_map %>%
#   write_csv('data/tranch1_relabel_map.csv')

tranch2_relabel_map = same_same_labelled$tranch2 %>%
  filter(match_level %in% c('mislabel', 'close')) %>%
  select(visit,subject,match_ratio) %>%
  dplyr::rename(same_same = match_ratio) %>%
  left_join(relabelled_samples_tranch2,
            by = c('subject', 'visit'),
            multiple='all') %>%
  mutate(dna_relabel = samples_greater_90,
         notes = "")

# tranch2_relabel_map %>%
#   write_csv('data/tranch2_relabel_map.csv')

