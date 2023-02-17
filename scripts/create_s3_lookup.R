library(tidyverse)
library(here)
library(RSQLite)

con = dbConnect(RSQLite::SQLite(), "/mnt/lts/personal/chasem/llfs/compile_database/data/pheno_data_202007.sqlite")

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

visit_split_bam_df = bam_df %>%
  group_by(visit) %>%
  group_split()

names(visit_split_bam_df) = unlist(map(visit_split_bam_df, ~unique(pull(.,visit))))

#' split a dataframe into chunk_size parts, so if there are 120 rows
#' and chunk_size is 60, you'd get two sets of 60 each. Last set might be
#' less than chunk_size. Names of each chunk in the returned list is the
#' iterated number of that chunk
create_chunks = function(subset_bam_df,chunk_size){

  n <- nrow(subset_bam_df)

  r <- rep(1:ceiling(n/chunk_size), each=chunk_size)[1:n]

  # for each of the subsets, copy the rows and rbind -- add the suffix
  # .bai to those copied rows. this will allow pulling both the
  # bam and the bai
   split(subset_bam_df,r) %>%
     map(~rbind(.,mutate(.,path = paste0(path,'.bai'))))

}

#' write chunks to remote mount. file prefix is what comes before the
#' chunk number, eg if visit_prefix is visit_2_chunk_ then the chunks
#' would get named visit_2_chunk_1 visit_2_chunk_2 ...
write_chunks = function(chunked_bam_df,
                        file_prefix,
                        output_dir ="/mnt/scratch/llfs_variant_calling/lookups"){

  map(names(chunked_bam_df),
    ~write_csv(select(chunked_bam_df[[.]],'path'),
               file.path(output_dir,
                         paste0(file_prefix,as.character(.))),
               col_names = FALSE))

}

visit_1_bam_split = create_chunks(visit_split_bam_df$visit1, 60)

visit_2_bam_split = create_chunks(visit_split_bam_df$visit2, 100)

write_chunks(visit_2_bam_split, 'visit_2_chunk_')



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

chunk='chunk_4'
bam_list = Sys.glob(paste0('/mnt/scratch/llfs_variant_calling/data/',chunk,'/*bam'))
visit_chunk_samplesheet = create_samplesheet_from_dir(bam_list, whatdatall, include_visit = TRUE)

write_csv(visit_chunk_samplesheet, paste0("/mnt/scratch/llfs_variant_calling/samplesheet/visit_1_",chunk,"_samplesheet.csv"))
#

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

create_sex_mislabels_comparison = function(row){

  row = as.data.frame(row)

  out = pedigree %>%
    filter(gpedid == row[['gpedid']],
           as.numeric(control) != as.numeric(row[['control']]))

  if(nrow(out) == 0){
    out = pedigree %>%
      filter(gpedid == row[['gpedid']])
  }

  out %>%
    dplyr::select(subject,gpedid,control) %>%
    dplyr::rename(dna_subject = subject, dna_control = control) %>%
    left_join(row %>%
                dplyr::select(subject,visit,gpedid) %>%
                dplyr::rename(rna_subject=subject),
              by='gpedid')
}

source("R/expand_grid_unique.R")

all_by_all_comparison = function(row, dna_subjects){
  rna_subject = row['subject']

  expand_grid_unique(rna_subject,dna_subjects) %>%
    as_tibble() %>%
    dplyr::rename(rna_subject = subject, dna_subject = z) %>%
    mutate(visit = row[['visit']])
}

comparison_lookup_maker = function(rna_dna_compare_df){

  rna_dna_compare_df %>%
    mutate(rna = file.path(
      '/scratch/mblab/chasem/llfs_rna_dna_compare_test',
      paste0(
        rna_subject,
        '_',
        visit,
        "_T1",
        ".haplotypecaller.filtered.gds"
      )
    )) %>%
    mutate(chr='1') %>%
    mutate(dna=file.path('/scratch/mblab/lisa.liao/human/staar/src/custom_scripts/agds',
                         paste0('LLFS.WGS.freeze5.chr',chr,'.gds'))) %>%
    dplyr::rename(rna_visit=visit) %>%
    mutate(rna_visit = str_remove(rna_visit,'visit')) %>%
    dplyr::select(rna_subject,rna_visit,dna_subject,rna,chr,dna)

}

sex_mislabel_tmp = map(seq(1,nrow(suspicious_sex_samples_set)),
                       ~create_sex_mislabels_comparison(suspicious_sex_samples_set[.,])) %>%
  do.call('rbind',.)

all_all_sex_mislabel_lookup = comparison_lookup_maker(sex_mislabel_tmp)

write_tsv(all_all_sex_mislabel_lookup,
          "/mnt/scratch/llfs_rna_dna_compare_test/sex_mislabel_all_comparison.tsv",
          col_names = FALSE)

dna_subjects =
  pedigree %>%
  filter(subject %in% bam_df$subject) %>%
  pull(subject)

sex_mislabel_tmp = map(seq(1,nrow(suspicious_sex_samples_set)),
                       ~all_by_all_comparison(suspicious_sex_samples_set[.,], dna_subjects)) %>%
  do.call('rbind',.)

sex_mislabel_lookup = sex_mislabel_tmp %>%
  mutate(rna = file.path(
    '/scratch/mblab/chasem/llfs_rna_dna_compare_test',
    paste0(
      rna_subject,
      '_',
      visit,
      "_T1",
      ".haplotypecaller.filtered.gds"
    )
  )) %>%
  mutate(count = rep(22,nrow(sex_mislabel_tmp))) %>%
  uncount(count) %>%
  group_by(rna_subject,dna_subject,visit) %>%
  mutate(chr=row_number(),
         dna=file.path('/scratch/mblab/lisa.liao/human/staar/src/custom_scripts/agds',
                       paste0('LLFS.WGS.freeze5.chr',row_number(),'.gds'))) %>%
  dplyr::rename(rna_visit=visit) %>%
  mutate(rna_visit = str_remove(rna_visit,'visit')) %>%
  dplyr::select(rna_subject,rna_visit,dna_subject,rna,chr,dna)

# write_tsv(sex_mislabel_lookup,
#           "/mnt/scratch/llfs_rna_dna_compare_test/sex_mislabel_compare_lookup.tsv",
#           col_names = FALSE)

