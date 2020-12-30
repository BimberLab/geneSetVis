testthat::test_that("Shiny app works", {
	shinytest::expect_pass(shinytest::testApp('../inst/app', compareImages = TRUE))
})
