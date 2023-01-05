
library(tidyverse)

create_samplesheet_from_dir = function(bam_list, whatdatall, include_visit=FALSE){
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
