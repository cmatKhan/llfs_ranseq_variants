library(tidyverse)
library(RSQLite)
library(SeqArray)
library(foreach)
library(here)

con = dbConnect(RSQLite::SQLite(), "~/projects/llfs_rnaseq_manifest/llfs_rnaseq_database.sqlite")

sample_df = tbl(con, 'correct_id_subject') %>%
  collect()

dna_chr21 = SeqArray::seqOpen(here('data/LLFS.WGS.freeze5.chr21.gds'))

wgs_sample_ids = seqGetData(dna_chr21, 'sample.id')

gds_files = list.files('/mnt/scratch/llfs_rna_dna_compare_test/new_samples_202307/data/',
                       '.*.gds',
                       recursive = TRUE,
                       full.names = TRUE)

rna_gds_df = tibble(rna_gds = gds_files,
                    batch_alias = basename(dirname(dirname(rna_gds))),
                    fastq_id = str_extract(basename(rna_gds), "^.*?(?=_T1)"),
                    visit = str_extract(fastq_id,'v\\d')) %>%
  replace_na(list(visit='0')) %>%
  mutate(fastq_id = str_remove(fastq_id, '_v3'),
         visit = str_extract(visit, "\\d"),
         rna_gds = str_replace(rna_gds, '/mnt/scratch/', '/scratch/mblab/chasem/')) %>%
  left_join(sample_df) %>%
  mutate(subject = ifelse(is.na(subject) & str_detect(fastq_id, 'pool'),
                          fastq_id, subject)) %>%
  mutate(wgs = subject %in% wgs_sample_ids) %>%
  filter(wgs)

all_by_all_full = foreach(
  i = seq(nrow(rna_gds_df)),
  .combine = 'bind_rows'
) %do% {
 row = rna_gds_df[i,]
 rna_gds = row[['rna_gds']]
 visit = row[['visit']]
 rna_subject = row[['subject']]
 rna_id = row[['fastq_id']]
 # note -- first, same same comparison, which is why rna_subject is the dna
 # subject. The RNA sample id is hte fastq_id
 tibble(rna_id = rna_id,
        visit = visit,
        dna_subject_id = rna_subject,
        #dna_subject_id = wgs_sample_ids,
        rna_gds = rna_gds,
        chr='1',
        dna_gds = '/ref/mblab/data/llfs/agds/LLFS.WGS.freeze5.chr1.gds')
}

all_by_all_full %>%
  write_tsv('/mnt/scratch/llfs_rna_dna_compare_test/lookups/all_by_all_post_20230519.txt',
            col_names = FALSE)

