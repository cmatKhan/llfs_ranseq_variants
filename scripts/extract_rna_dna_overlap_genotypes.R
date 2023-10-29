#!/usr/bin/env Rscript

# Dependencies ----------------------------------------------------------------
library(optparse)
library(SeqArray)
library(SeqVarTools)
library(foreach)
library(tidyverse)

# local testing ----------------------------------------------------------------
TEST = FALSE

# helper functions -------------------------------------------------------------
#
#' set chromosome and sample filter on gds obj
#'
#' @param gds_obj an open gds obj
#' @param chr chromosome number (no prefix)
#' @param sample_id name of the sample
set_initial_gds_filters = function(gds_obj,chr,sample_id){
  if(!as.character(sample_id) %in% seqGetData(gds_obj, "sample.id")){
    stop(sprintf("sample id %s does not exist in %s",sample_id,gds_obj$filename))
  } else if(!as.character(chr) %in% names(table(seqGetData(gds_obj, "chromosome")))){
    stop(sprintf("Chromosome %s not in %s",as.character(chr),gds_obj$filename))
  }
  seqSetFilterChrom(gds_obj,as.numeric(chr))
  seqSetFilter(gds_obj,sample.id=as.character(sample_id))
  set_variant_filter(gds_obj)

  length(nAlleles(gds_obj))
}

#' set variant level filters, currently such that only variants with 1
#' alternate, and only single nucleotide variant locations, are considered
#' @param gds_obj an open gds object
set_variant_filter = function(gds_obj){

  # select only bi-allelic sites
  sel_biallelic = which(seqNumAllele(gds_obj) == 2)

  # select only snv sites
  snv_variants = which(isSNV(gds_obj))

  sel_variants = intersect(sel_biallelic,snv_variants)
  seqSetFilter(gds_obj, variant.sel = sel_variants, action='intersect')
}

#' find variant position overlap between rna and dna, set filter accordingly
#' @note this could be generalized to find overlap between any number of gds
#' objects in the list
#' @param gds_list right now, explicitly written for a list with two slots named
#' rna and dna
set_overlap_filter = function(gds_list){
  # extract granges
  rna_granges = granges(gds_list$rna)
  dna_granges = granges(gds_list$dna)
  # find overlap
  var_overlap = IRanges::subsetByOverlaps(rna_granges,dna_granges,type='equal')

  # set filters
  seqSetFilterPos(
    gds_list$rna,
    as.numeric(opt$chr),
    pos=var_overlap@ranges@start)

  seqSetFilterPos(
    gds_list$dna,
    as.numeric(opt$chr),
    pos=var_overlap@ranges@start)

  length(nAlleles(gds_list$rna))
}

#' extract genotype matrix, tidy and return as tibble
#'
#' @param gds_obj open gds object
#' @param source  eg 'rna' or 'dna'. This column stores the genotype
#' @param sample eg 1085, the subject/sample id
tidy_geno_mat = function(gds_obj,source,sample,chr){
  source_sample_colname = paste(source,'id',sep="_")
  # see here for diff btwn info/DP and format/DP
  # https://gatk.broadinstitute.org/hc/en-us/articles/360035532252-Allele-Depth-AD-is-lower-than-expected
  total_depth_colname = paste(source,'total_depth',sep='_')
  total_depth = seqGetData(gds_obj,'annotation/format/DP')[1,]
  ad = seqGetData(gds_obj,'annotation/format/AD')
  alt_depth_colname = paste(source,'alt_depth',sep="_")
  alt_depth = ad$data[1,seq(2,length(ad$length)*2,2)]
  # see here for GQ vs PL
  # https://gatk.broadinstitute.org/hc/en-us/articles/360035890451-Calculation-of-PL-and-GQ-by-HaplotypeCaller-and-GenotypeGVCFs
  genotype_qual_colname = paste(source,'gq',sep="_")
  genotype_qual = seqGetData(gds_obj,'annotation/format/GQ')[1,]

  getGenotype(gds_obj) %>%
    as_tibble() %>%
    pivot_longer(everything(),values_to=source,names_to = 'var_id') %>%
    mutate(!!rlang::sym(source) := str_replace(!!rlang::sym(source),"\\|","/"),
           !!rlang::sym(source_sample_colname) := sample) %>%
    left_join(
      granges(gds_obj) %>%
        as_tibble() %>%
        mutate(var_id = names(granges(gds_obj))) %>%
        dplyr::select(var_id,start) %>%
        dplyr::rename(pos=start)) %>%
    mutate(chr=chr) %>%
    dplyr::select(chr,pos,all_of(source_sample_colname),all_of(source)) %>%
    mutate(!!rlang::sym(total_depth_colname) := total_depth,
           !!rlang::sym(alt_depth_colname) := alt_depth,
           !!rlang::sym(genotype_qual_colname) := genotype_qual)
}

# cmd line argument parser -----------------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c("-v", "--verbose"), action="store_true",
                     default=TRUE, help="Print extra output [default]")

parser <- add_option(parser, c("-q", "--quietly"), action="store_false",
                     dest="verbose", help="Print little output")

parser <- add_option(parser, c("--chr"), type="character",
                     default="",
                     help="chromomsome on which to extract overlap data",
                     metavar="'chr1'")

parser <- add_option(parser, c("--rna"), type="character",
                     default="",
                     help="path to rna gds file",
                     metavar="'/path/to/rna.gds'")

parser <- add_option(parser, c("--dna"), type="character",
                     default="",
                     help="path to dna gds file",
                     metavar="'/path/to/dna.gds'")

parser <- add_option(parser, c("--write_full_comparison"),
                     action = "store_true",
                     help="set this to output both the full sample to sample
                      comparison and the similarity summary",
                     default=FALSE)

parser <- add_option(parser, c('--rna_sample'), type="character",
                    default="",
                    help="sample id used to extract data from the rna gds files",
                    metavar ="'sample.id'")

parser <- add_option(parser, c('--rna_visit'), type="character",
                    default="1",
                    help="visit number of rna sample",
                    metavar ="'1'")

parser <- add_option(parser, c('--dna_sample'), type="character",
                    default="",
                    help="sample id used to extract data from the dna gds files",
                    metavar="'sample.id'")

parser <- add_option(parser, c("--output_prefix"), type="character",
                     default=".",
                     help="output location",
                     metavar="'.'")

# main -------------------------------------------------------------------------

opt = if(TEST){
  parse_args(parser,
           args = c(
             "--quietly",
             "--chr=1",
             "--rna_sample=20251",
             "--dna_sample=20251",
             "--dna=/mnt/ref/data/llfs/geno_chip/llfs_gwas.chr1.gds",
             "--rna=/mnt/scratch/llfs_rna_dna_compare_test/tmp_vcf_gds/20251_visit2_T1.haplotypecaller.filtered.gds"))
} else{
  parse_args(parser)
}

# check cmd line arguments
x = foreach(
  i = names(opt)
) %do% {
  input_value=opt[[i]]
  if(input_value==''){
    stop(sprintf("ARGUMENT --%s IS REQUIRED",i))
  }else if(i %in% c('dna','rna')){
    if(!file.exists(input_value)){
      stop(sprintf("FILE %s DOES NOT EXIST",input_value))
    }
  }
}

if(opt$verbose){
  message(sprintf("rna_sample: %s; dna_sample: %s; chr: %s",
                opt$rna_sample,opt$dna_sample,opt$chr))
}


# uncomment this to profile the process
# p <- profvis({
# memory profiling result for chr21 shows max usage at 350MB
# time at 36 seconds

# open gds files
gds_list = list(
  rna = seqOpen(opt$rna,allow.duplicate = TRUE),
  dna = seqOpen(opt$dna,allow.duplicate = TRUE)
)

init_rna_var_n = set_initial_gds_filters(
  gds_list$rna,
  opt$chr,
  opt$rna_sample)

init_dna_var_n = set_initial_gds_filters(
  gds_list$dna,
  opt$chr,
  opt$dna_sample)

# apply position filter based on overlapping variants
overlap_var_n = set_overlap_filter(gds_list)

# extract data
compare_df = tidy_geno_mat(gds_list$rna,'rna',opt$rna_sample,opt$chr) %>%
  left_join(tidy_geno_mat(gds_list$dna,'dna',opt$dna_sample,opt$chr),by=c('chr','pos')) %>%
  filter(!is.na(rna),!is.na(dna),rna_total_depth <= 300) %>%
  mutate(homo_match = ((rna=='1/1' & dna=='1/1') | (rna=='0/0' & dna=='0/0')),
         hetero_match = (rna=='0/1'& dna=='0/1'),
         homo_expr_candidate = ((rna=='1/1' & dna=='0/1') | (rna=='0/0' & dna=='0/1')))

comparison_summary_unfltr = compare_df %>%
  group_by(chr,homo_match,hetero_match,homo_expr_candidate) %>%
  tally() %>%
  mutate(dna_sample = opt$dna_sample,rna_sample = opt$rna_sample, depth_filter=0,gq_filter=0)

depth_fltr=10
gq_fltr=80
comparison_summary_fltr = compare_df %>%
  filter(dna_total_depth >= depth_fltr,
         rna_total_depth >= depth_fltr,
         rna_gq >= gq_fltr,
         dna_gq >= gq_fltr) %>%
  group_by(homo_match,hetero_match,homo_expr_candidate) %>%
  tally() %>%
  mutate(dna_sample = opt$dna_sample,
         rna_sample = opt$rna_sample,
         depth_filter=depth_fltr,
         gq_filter = gq_fltr)

unfilter_match_n = comparison_summary_unfltr %>%
  filter(homo_match==TRUE | hetero_match==TRUE) %>%
  pull(n) %>%
  sum()

unfiltered_match_ratio = unfilter_match_n / sum(comparison_summary_unfltr$n) %>%
  round(.,digits=4)

homo_expr_cand_unfltr = comparison_summary_unfltr %>%
  filter(homo_expr_candidate == TRUE) %>%
  pull(n) %>%
  ifelse(identical(.,integer(0)),0,.)

filtered_match_n = comparison_summary_fltr %>%
  filter(homo_match==TRUE | hetero_match==TRUE) %>%
  pull(n) %>%
  sum()

filtered_match_ratio = filtered_match_n / sum(comparison_summary_fltr$n) %>%
  round(.,digits=4)

homo_expr_cand_fltr = comparison_summary_fltr %>%
  filter(homo_expr_candidate == TRUE) %>%
  pull(n) %>%
  ifelse(identical(.,integer(0)),0,.)

sample_similarity_metric = tibble(
  chr = opt$chr,
  rna_sample = opt$rna_sample,
  rna_visit = opt$rna_visit,
  dna_sample = opt$dna_sample,
  rna_var_n = init_rna_var_n,
  dna_var_n = init_dna_var_n,
  overlap_var_n = overlap_var_n,
  overlap_unfltr = sum(comparison_summary_unfltr$n),
  overlap_fltr = sum(comparison_summary_fltr$n),
  n_match_unfltr = unfilter_match_n,
  n_match_fltr = filtered_match_n,
  homo_expr_cand_fltr = homo_expr_cand_fltr
)


output_dir  = paste('rna',opt$rna_sample,sep="_")
system_output_dir = file.path(opt$output_prefix,output_dir)
dir.create(system_output_dir,recursive=TRUE,showWarnings=TRUE)

file_name = paste('visit', opt$rna_visit, 'chr',opt$chr, opt$dna_sample, sep = "_")

# write full comparison df if -f is set
if(opt$write_full_comparison){
  write_csv(compare_df,
            file.path(system_output_dir,
                      paste0(file_name,"_comparison.csv")))
}
# write summarized match metric
write_csv(sample_similarity_metric,
          file.path(system_output_dir,
                    paste0(file_name,"_match_metrics.csv")))

# this is part of profvis -- uncomment this for resource profiling during testing
# })
# print(p)

