library(tidyverse)
library(here)
library(RSQLite)

con = dbConnect(RSQLite::SQLite(), "~/projects/llfs/compile_database/data/pheno_data_202007.sqlite")

whatdatall = tbl(con,"whatdatall") %>%
  collect()

pedigree = read_csv("data/triplet_visit2_version3.csv")

df = read_tsv("/mnt/scratch/s3_info/llfs_20221207_bam_only.txt",col_names = 'path')

suspicious_sex_samples = read_csv("data/suspicious_sex_samples.csv") %>%
  left_join(whatdatall %>%
              dplyr::select(subject,id) %>%
              distinct(subject,.keep_all = TRUE),by='subject') %>%
  mutate(visit = str_remove(str_extract(visit,'visit\\d|visit_\\d|pool'),"_"))

bam_df = df %>%
  mutate(visit = str_remove(str_extract(path,'visit\\d|visit_\\d|pool'),"_"),
         id = str_extract(basename(path), "\\d+")) %>%
  left_join(
    select(whatdatall,id,subject) %>%
      distinct(id,.keep_all = TRUE) %>%
      filter(!is.na(subject))) %>%
  left_join(pedigree)

suspicious_sex_samples_set = suspicious_sex_samples %>%
  left_join(bam_df,by=c('id','subject','visit','sex'))

bam_df %>%
  group_by(visit,gpedid) %>%
  tally()

sample_from_df = bam_df %>%
  filter(visit == 'visit1', gen %in% c(1,2,3)) %>%
  group_by(gpedid) %>%
  filter(n() >=40)

related_samples = sample_from_df %>%
  filter(relative == 1)

control_samples = sample_from_df %>%
  filter(control == 1)

stopifnot(setdiff(related_samples$gpedid, control_samples$gpedid) ==
            setdiff(control_samples$gpedid,related_samples$gpedid))

test_set = rbind(slice_sample(control_samples,n=3),slice_sample(related_samples,n=3))

test_set %>%
  ungroup() %>%
  select(path)
# %>%
#   write_tsv("/mnt/scratch/llfs_variant_calling/test_set_s3.txt", col_names = FALSE)

create_samplesheet_from_dir = function(bam_list, include_visit=FALSE){
  df = tibble(bam=bam_list) %>%
    mutate(bam = str_replace(bam,"/mnt/scratch","/scratch/mblab/chasem"),
           visit = str_remove(str_extract(bam,'visit\\d|visit_\\d|pool'),"_"),
           id = str_extract(basename(bam), "\\d+")) %>%
    left_join(select(whatdatall,id,subject) %>%
                distinct(id,.keep_all = TRUE) %>%
                filter(!is.na(subject))) %>%
    mutate(fastq_1="",fastq_2="",strandedness="") %>%
    dplyr::rename(sample = subject)

  if(include_visit){
    df = df %>%
      mutate(sample = paste(sample,visit,sep='_'))
  }
    dplyr::select(df,sample,fastq_1,fastq_2,strandedness,bam)
}

bam_list = Sys.glob('/mnt/scratch/llfs_variant_calling/data/sex_mislabels/*bam')
sex_mislabel_samplesheet = create_samplesheet_from_dir(bam_list, include_visit = TRUE)

write_csv(sex_mislabel_samplesheet, "/mnt/scratch/llfs_variant_calling/sex_mislabel_samplesheet.csv")

# samplesheet_df = tibble(bam=Sys.glob("/mnt/scratch/llfs_variant_calling/data/llfs_test_set/*bam")) %>%
#   mutate(bam = str_replace(bam,"/mnt/scratch","/scratch/mblab/chasem"),
#          visit = str_remove(str_extract(bam,'visit\\d|visit_\\d|pool'),"_"),
#          id = str_extract(basename(bam), "\\d+")) %>%
#   left_join(select(whatdatall,id,subject) %>%
#               distinct(id,.keep_all = TRUE) %>%
#               filter(!is.na(subject))) %>%
#   mutate(fastq_1="",fastq_2="",strandedness="") %>%
#   dplyr::rename(sample = subject) %>%
#   dplyr::select(sample,fastq_1,fastq_2,strandedness,bam)



# write_csv(samplesheet_df,"/mnt/scratch/llfs_variant_calling/llfs_test_samplesheet.csv")

samplesheet_df = read_
