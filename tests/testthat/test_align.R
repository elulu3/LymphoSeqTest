context("Perform multiple sequence alignment")
library(LymphoSeq2)

test_that("Align all sequences in all sample within edit distance of 15", {
  base::set.seed(12357)
  stable <- LymphoSeq2::readImmunoSeq("test_data/")
  ntable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction")
  nalign <- LymphoSeq2::alignSeq(ntable)
  nseq <- base::length(nalign@unmasked)
  known_consensus <- "---------?T??AG?CC?C?GAGC??G??GACTC?GCC?TGTAT?TCTGTGCCAGCAGC???G???????????????????????---------??????TT?GG????---"
  test_consensus <- msa::msaConsensusSequence(nalign)
  expect_equal(nseq, 265)
})

test_that("Align all sequences in one sample within edit distance of 15", {
  ttable <- LymphoSeq2::readImmunoSeq("test_data/")
  tntable <- LymphoSeq2::productiveSeq(ttable, aggregate = "junction")
  talign <- LymphoSeq2::alignSeq(tntable, repertoire_ids = "015V06013979_CFAR")
  tseq <- base::length(talign@unmasked)
  tconsensus <- msa::msaConsensusSequence(talign)
  tname <- base::unique(base::names(talign@unmasked))
  ktable <- LymphoSeq2::readImmunoSeq("test_data/015V06013979_CFAR.tsv")
  kntable <- LymphoSeq2::productiveSeq(ktable, aggregate = "junction")
  kalign <- LymphoSeq2::alignSeq(kntable)
  kseq <- base::length(kalign@unmasked)
  kconsensus <- msa::msaConsensusSequence(kalign)
  kname <- "015V06013979_CFAR"
  expect_equal(tconsensus, kconsensus)
  expect_equal(tseq, kseq)
  expect_equal(tname, kname)
})

test_that("Align single sequence in all samples within edit distance of 15", {
  base::set.seed(12357)
  ttable <- LymphoSeq2::readImmunoSeq("test_data/")
  tntable <- LymphoSeq2::productiveSeq(ttable, aggregate = "junction")
  talign <- LymphoSeq2::alignSeq(tntable, sequence_list = c("AACGCCTTGGAGCTGGACGACTCGGCCATGTATCTCTGTGCCAGCAGCTTGGCGGGGGGGCCGTGGGAGACCCAGTACTTCGGGCCA"))
  tseq <- base::length(talign@unmasked)
  tconsensus <- msa::msaConsensusSequence(talign)
  tname <- length(base::unique(base::names(talign@unmasked)))
  kseq <- 6
  kconsensus <- "CA?CCC???GAGCTGGAGGACTCGGCC?TGTATCTCTGTGCCAGCAGCTTGGCGGGGGGGCCGTGGGAGACCCAGTACTTCGGGCCA"
  kname <- 2
  expect_equal(tconsensus, kconsensus)
  expect_equal(tseq, kseq)
  expect_equal(tname, kname)
})
