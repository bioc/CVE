library(testthat)
library(CVE)

test_that("oncotator example file is in the right format",
          {
            expect_that(oncotator_example, is.data.frame)
          })


