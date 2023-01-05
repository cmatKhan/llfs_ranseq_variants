library(SeqArray)
library(SeqVarTools)
library(tidyverse)
library(parallel)
library(here)

file = seqOpen(here("data/LLFS.WGS.freeze5.chr21.gds"))

depth_summary_by_sample = function(sample){
  seqSetFilter(file,sample.id=sample)
  broom::tidy(summary(seqGetData(file,'annotation/format/DP')))
  seqResetFilter(file)
}

AF <- seqAlleleFreq(file)
seqSetFilter(file,variant.sel=(AF<=0.5))

sel_by_qual = which(qual(file) > 20)
seqSetFilter(file,sel_by_qual)

# select only bi-allelic sites
sel_biallelic = which(seqNumAllele(file) == 2)
seqSetFilter(file,sel_biallelic)

