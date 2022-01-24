context("Find n top sequences with highest frequency")
library(LymphoSeqTest)


test_that("Find top sequence in all samples", {
  stable <- LymphoSeqTest::readImmunoSeq("test_data/015V06013979_CFAR.tsv")
  atable <- LymphoSeqTest::productiveSeq(stable, aggregate = "junction_aa")
  ttable <- LymphoSeqTest::topSeqs(atable, top = 1) %>%
            dplyr::pull(duplicate_frequency) %>%
            base::max()
  ctable <- atable %>% 
            dplyr::pull(duplicate_frequency) %>%
            base::max()
  expect_equal(ctable, ttable)
})

test_that("Find the top ten sequences in one sample", {
  stable <- LymphoSeqTest::readImmunoSeq("test_data/015V06013979_CFAR.tsv")
  atable <- LymphoSeqTest::productiveSeq(stable, aggregate = "junction_aa")
  ttable <- LymphoSeqTest::topSeqs(atable, top = 10) %>%
            dplyr::pull(duplicate_frequency) 
  ctable <- atable %>% 
            dplyr::pull(duplicate_frequency) %>%
            base::unique() %>%
            base::sort(decreasing = TRUE)
  expect_equal(ctable[1:10], ttable)
})


test_that("Top sequences across samples is correctly identified", {
  stable <- LymphoSeqTest::readImmunoSeq("test_data/")
  atable <- LymphoSeqTest::productiveSeq(stable, aggregate = "junction_aa")
  ttable <- LymphoSeqTest::topSeqs(atable, top = 1) 
  
  ctable <- atable %>% 
            dplyr::group_by(repertoire_id) %>%
            dplyr::group_split(.keep = TRUE) %>%
            purrr::map(LymphoSeq2::topSeqs) %>%
            dplyr::bind_rows()
  expect_true(base::all.equal(ctable, ttable))
})


