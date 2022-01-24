context("Read GLIPH Files")
library(LymphoSeqTest)

test_that("Reads a GLIPH files correctly", {
  gtable <- LymphoSeqTest::readGliph("test_data/gliph")
  nsample <- gtable %>%
             dplyr::pull(repertoire_id) %>%
             base::unique() %>% 
             base::length()
  nspec <- gtable %>%
           dplyr::pull(spec_group) %>%
           base::unique() %>%
           base::length()
  expect_equal(nsample, 2)
  expect_equal(nspec, 20)
})
