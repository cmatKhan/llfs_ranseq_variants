library(GenomicFeatures)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

variants = list(
  dna_variants = readRDS("data/dna_granges_20251_genoarray.rds"),
  rna_variants = readRDS("data/rna_granges_20251.rds")
)

features = list(
  introns = intronicParts(txdb),
  exons = exonicParts(txdb),
  # https://support.bioconductor.org/p/66003/#66046
  intergenic = gaps(genes(txdb))
)
# make the seqnames the same as the annotations
variants$dna_variants@seqinfo@seqnames = 'chr1'
variants$rna_variants@seqinfo@seqnames = 'chr1'

get_feature_overlaps = function(var_ranges){
  map(features, ~length(findOverlaps(var_ranges,.)))
}

map(variants, get_feature_overlaps) %>%
  as_tibble() %>%
  mutate(dna_variants = as.numeric(dna_variants),
         rna_variants = as.numeric(rna_variants),
         feature = names(features)) %>%
  pivot_longer(-feature,names_to='source', values_to='count') %>%
  ggplot(aes(feature, count, fill=source, group=source)) +
  geom_bar(position='dodge', stat='identity') +
  theme_bw()+
  theme(text = element_text(size=20))

