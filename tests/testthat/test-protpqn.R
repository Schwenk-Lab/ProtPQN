test_that("correct output format", {
  protein_df <- replicate(50, rnorm(100, 10, 1))
  expect_equal(class(protpqn(protein_df)), class(protein_df))
  expect_equal(dim(protpqn(protein_df)), dim(protein_df))

  protein_df <- as.data.frame(protein_df)
  expect_equal(class(protpqn(protein_df)), class(protein_df))
  expect_equal(dim(protpqn(protein_df)), dim(protein_df))
})

test_that("error for inf", {
  protein_df <- replicate(50, rnorm(100, 10, 1))
  # Introduce Inf anywhere
  protein_df[sample(nrow(protein_df), 1), sample(ncol(protein_df), 1)] <- Inf

  expect_error()
})
