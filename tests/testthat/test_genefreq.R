context("Get VDJ gene frequency")
library(LymphoSeqTest)


test_that("Gene frequency charcterized by gene name sums to 1", {
  stable <- LymphoSeqTest::readImmunoSeq("test_data")
  ntable <- LymphoSeqTest::productiveSeq(stable, aggregate = "junction")
  gtable <- LymphoSeqTest::geneFreq(ntable) %>% 
            dplyr::group_by(repertoire_id, gene_type) %>% 
            dplyr::summarise(freq_tot = sum(gene_frequency)) %>%
            dplyr::pull(freq_tot)
  expect_true(base::all(gtable == 1))
})

test_that("Gene frequency charcterized by gene name sums to 1", {
  stable <- LymphoSeqTest::readImmunoSeq("test_data")
  ntable <- LymphoSeqTest::productiveSeq(stable, aggregate = "junction")
  gtable <- LymphoSeqTest::geneFreq(ntable, family = TRUE) %>% 
            dplyr::group_by(repertoire_id, gene_type) %>% 
            dplyr::summarise(freq_tot = sum(gene_frequency)) %>%
            dplyr::pull(freq_tot)
  expect_true(base::all(gtable == 1))
})
    
