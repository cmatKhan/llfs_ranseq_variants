library(tidyverse)
library(here)
library(RSQLite)

source(here("R/create_samplesheet_from_dir.R"))

con = dbConnect(RSQLite::SQLite(), "/mnt/lts/personal/chasem/llfs/compile_database/data/pheno_data_202007.sqlite")

whatdatall = tbl(con,"whatdatall") %>%
  collect()

pedigree = read_csv("data/triplet_visit2_version3.csv")

chunk='chunk_11'
bam_list = Sys.glob(paste0('/mnt/scratch/llfs_variant_calling/data/',chunk,'/*bam'))
visit_chunk_samplesheet = create_samplesheet_from_dir(
  bam_list,
  whatdatall,
  include_visit = TRUE) %>%
  mutate(sample = ifelse(
    str_detect(sample,"^NA"),
    paste0(str_extract(basename(bam), "^\\d+"), "_visit2"), sample))

write_csv(visit_chunk_samplesheet, paste0("/mnt/scratch/llfs_variant_calling/samplesheet/visit_2_",chunk,"_samplesheet.csv"))

dbDisconnect(con)

