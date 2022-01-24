context("Calculate clonal realtedness of samples")
library(LymphoSeqTest)

test_that("Calculate clonal relatedness of all sequences", {
  stable <- LymphoSeqTest::readImmunoSeq("test_data/")
  ttable <- LymphoSeqTest::clonalRelatedness(stable)
  ctable <- tibble::tibble(repertoire_id = c("015V06013979_CFAR", "015V12001549_CFAR", "015V12001685_CFAR_R", "015V12003105_CFAR"),
                           clonalRelatedness = c(0.02247191, 0.02000000, 0.02000000, 0.02000000))
  expect_true(base::all.equal(ttable, ctable))
})

test_that("Calculate clonal relatedness of all productive nucleotide sequences", {
  stable <- LymphoSeqTest::readImmunoSeq("test_data/")
  ntable <- LymphoSeqTest::productiveSeq(stable, aggregate = "junction")
  ttable <- LymphoSeqTest::clonalRelatedness(ntable) %>%
            dplyr::mutate(clonalRelatedness = base::round(clonalRelatedness, 3))
  ctable <- tibble::tibble(repertoire_id = c("015V06013979_CFAR", "015V12001549_CFAR", "015V12001685_CFAR_R", "015V12003105_CFAR"),
    clonalRelatedness = c(0.022, 0.025, 0.026, 0.027))
  expect_true(base::all.equal(ttable, ctable))
})
