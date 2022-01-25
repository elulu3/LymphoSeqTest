context("Count k-mers in nucleotide sequence")
library(LymphoSeqTest)

test_that("Number of k-mers are counted correctly", {
    junction <- "ACCTAGGT"
    study_table <- tibble(junction)
    ktable <- LymphoSeqTest::countKmer(study_table = study_table, k = 5)
    num_rows <- base::nrow(ktable)
    expect_equal(num_rows, 4)
})