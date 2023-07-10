library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(SeqVarTools)
library(tidyverse)
library(RSQLite)
library(readxl)

WRITE=FALSE
ONLY_NEW=TRUE

# generate list of WGS sample ids
llfs_gds <- seqOpen("data/LLFS.WGS.freeze5.chr21.gds")
wgs_sample_ids <- seqGetData(llfs_gds, "sample.id")
seqClose(llfs_gds)
# generate list of geno array sample ids
# NOTE: it turns out geno arr vcf are hg19
geno_arr_gds <- seqOpen("/mnt/ref/data/llfs/geno_chip/original_hg19/llfs_gwas.chr22.gds")
geno_arr_sample_id <- seqGetData(geno_arr_gds, "sample.id")
seqClose(geno_arr_gds)

# generate the sample to id map
subj_id_map <- read_csv("data/whatdatall.csv") %>%
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

# parse the phantom ID samples
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
  left_join(subj_id_map, by = c("id" = "whatdatall_id")) %>%
  mutate(phan_subj = ifelse(is.na(phan_subj), whatdatall_subj, phan_subj))

# parse the legal fails
legal_fails <- read_excel(
  file.path(
    "/mnt/lts/personal/chasem/llfs",
    "data_processing/data/LLFS_DKsamples_NOT_RNASeq_20220523.xlsx"
  )
)

s3_vcf_path = "data/all_filtered_vcf_20230313.txt"
s3_vcf_all <- read_csv(s3_vcf_path, col_names = "vcf") %>%
  filter(str_detect(vcf, "tbi$", negate = TRUE)) %>%
  mutate(
    subject = str_remove(
      basename(vcf),
      "_visit\\d_T1.haplotypecaller.filtered.vcf.gz.*"
    ),
    visit = str_extract(vcf, "visit_\\d|visit\\d"),
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

if(ONLY_NEW){
  # filter down to only the newly sequenced samples
  s3_vcf_all = s3_vcf_all %>%
    filter(str_detect(vcf, "s6418|s6419|s6420"))
}

s3_vcf_all_not_in_dna <- filter(s3_vcf_all, !in_wgs_dna, !in_genoarr_dna) %>%
  filter(str_detect(subject, 'pool', negate = TRUE))

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

if(WRITE){
  full_same_same_comp %>%
    # filter(subject %in% c(13757, 14752)) %>%
  write_tsv(
    "/mnt/scratch/llfs_rna_dna_compare_test/lookups/full_same_same_comp_makeup.txt",
    col_names = FALSE)
}

# full_same_same_comp %>%
#   filter(subject %in% c(13757, 14752)) %>%
#   mutate(gds = str_replace(vcf,'vcf.gz','gds')) %>%
#   select(subject, visit, gds) %>%
#   write_tsv('/mnt/scratch/llfs_rna_dna_compare_test/new_samples_makeup/lookups/rna_samples.txt',
#             col_names=FALSE)

if(WRITE){
  s3_vcf_all_not_in_dna %>%
    select(-c(vcf, id, phan_subj, missing_subject, whatdatall_subj)) %>%
    write_csv("data/rna_subjects_not_in_dna_20230217.csv")
}
