
#setwd('app')

# to create expected results; can only save under app/tests so move to scripts to ../../test after
#shinytest::recordTest("inst/app", loadTimeout=100000, seed = 1234)

testthat::test_that("Application works", {
  # Use compareImages=FALSE because the expected image screenshots were created
  # on a Mac, and they will differ from screenshots taken on the CI platform,
  # which runs on Linux.
  #shinytest::expect_pass(snapshotPreprocessOutput(shinytest::testApp('inst/app', compareImages = FALSE, quiet = F)))
  #shinytest::expect_pass(shinytest::testApp('inst/app', compareImages = FALSE, quiet = F))
  shinytest::expect_pass(shinytest::testApp('./inst/app', compareImages = FALSE))
})
##ERROR:Unknown whether app is a regular Shiny app or .Rmd: ./inst/app
