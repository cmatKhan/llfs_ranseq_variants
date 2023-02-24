#!/usr/bin/env Rscript

library(futile.logger)
library(optparse)
library(furrr)
library(future)
library(tidyverse)

# GLOBALS ---------------------------------------------------------------------

TEST = FALSE

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
parser <- add_option(parser, c("--cpus"), type="integer",
                     default=10,
                     help="number of CPUs",
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
               "--cpus=10",
               "--dir=/mnt/scratch/test_dir"))
} else{
  parse_args(parser)
}

futile.logger::flog.info(sprintf('setting up parallel env with %s cpus', opt$cpus))

future::plan(
  future::multisession,
  workers=as.integer(opt$cpus)
)


gather_data_from_all_by_all = function(file_list){


  furrr::future_map(file_list,read_csv) %>%
    do.call('rbind',.)
}

metrics_file_list = list.files(
  opt$dir,
  "*_match_metrics.csv",
  recursive = TRUE,
  full.names = TRUE)

futile.logger::flog.info(sprintf('reading in data from %s files', length(metrics_file_list)))
df = gather_data_from_all_by_all(metrics_file_list)

write_csv(df,opt$output)



