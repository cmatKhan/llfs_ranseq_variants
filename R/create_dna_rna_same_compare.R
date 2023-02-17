library(tidyverse)

x = read_tsv("/mnt/scratch/llfs_rna_dna_compare_test/variant_lookup.tsv",
             col_names = FALSE)


parsed_lookup = filter(x,str_detect(X4,"gvcf|gds",negate = TRUE)) %>%
  filter(str_detect(X4, "filtered")) %>%
  select(X4) %>%
  mutate(subject = str_extract(basename(X4),"^\\d+")) %>%
  dplyr::rename(s3_path=X4) %>%
  filter(str_detect(s3_path,".tbi$",negate = TRUE)) %>%
  mutate(index = paste0(s3_path,".tbi")) %>%
  select(s3_path,index,subject)

compare_results_df = readRDS("data/intermediate_all_same_compare_20230124.rds")

parsed_lookup %>%
  filter(!subject %in% compare_results_df$rna_sample) %>%
  write_tsv("/mnt/scratch/llfs_rna_dna_compare_test/lookups/finish_visit1_rna_dna.txt",
            col_names = FALSE)
