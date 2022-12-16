#!/usr/bin/env Rscript

library(optparse)
library(SeqArray)
library(tidyverse)
library(foreach)


# arg parser -------------------------------------------------------------------
parser <- OptionParser()
parser <- add_option(parser, c("-v", "--verbose"), action="store_true",
                     default=TRUE, help="Print extra output [default]")

parser <- add_option(parser, c("-q", "--quietly"), action="store_false",
                     dest="verbose", help="Print little output")


parser <- add_option(parser, c("--vcf"), type="character",
                     default="",
                     help="path to vcf file",
                     metavar="'/path/to/file.vcf'")

parser <- add_option(parser, c("--output_prefix"), type="character",
                     default=".",
                     help="output location",
                     metavar="'.'")

# main -------------------------------------------------------------------------

opt = parse_args(parser)

# check cmd line arguments
x = foreach(
  i = names(opt)
) %do% {
  input_value=opt[[i]]
  if(input_value==''){
    stop(sprintf("ARGUMENT --%s IS REQUIRED",i))
  }else if(i %in% c('vcf')){
    if(!file.exists(input_value)){
      stop(sprintf("FILE %s DOES NOT EXIST",input_value))
    }
  }
}

gds_filename = str_replace_all(basename(opt$vcf),'.vcf.gz|.vcf','.gds')

seqVCF2GDS(opt$vcf,gds_filename)
