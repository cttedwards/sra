test_that("compile works", {
  model <- sra()
  expect_s4_class(model, "stanmodel")
})
