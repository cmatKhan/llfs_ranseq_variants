#!/usr/bin/env Rscript

library(optparse)
library(tidyverse)
library(furrr)
library(future)

# GLOBALS ---------------------------------------------------------------------

TEST = TRUE

# cmd line argument parser ----------------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c("-v", "--verbose"), action="store_true",
                     default=TRUE, help="Print extra output [default]")
parser <- add_option(parser, c("-q", "--quietly"), action="store_false",
                     dest="verbose", help="Print little output")
parser <- add_option(parser, c("--dir"), type="character",
                     default=".",
                     help="directory to search",
                     metavar ="'directory'")
parser <- add_option(parser, c("--output"), type="character",
                     default="compiled_metrics.csv",
                     help="name (may include path, otherwise PWD) of file to output. output format is csv",
                     metavar ="'output'")


# main -------------------------------------------------------------------------

opt = if(TEST){
  parse_args(parser,
             args = c(
               "--quietly",
               "--dir=/mnt/scratch/test_dir"))
} else{
  parse_args(parser)
}

future::plan(
  future::multisession,
  workers=10
)

gather_data_from_all_by_all = function(){
  compare_results = list.files(
    opt$dir,
    "*_match_metrics.csv",
    recursive = TRUE,
    full.names = TRUE)

  furrr::future_map(compare_results,read_csv) %>%
    do.call('rbind',.)
}

df = gather_data_from_all_by_all()

write_csv(df,opt$output)



