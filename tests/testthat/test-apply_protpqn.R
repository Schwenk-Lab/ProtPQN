test_that("correct output format", {
  # Wide format test data
  protein_data <- replicate(50, rnorm(100, 10, 1)) |> as.data.frame() |> log2()
  proteins <- paste0("p", 1:50)
  colnames(protein_data) <- proteins
  protein_data$sample_id <- paste0("s", 1:100)
  protein_data <- protein_data[, c(ncol(protein_data), 1:(ncol(protein_data)-1))]

  normalized_data <- apply_protpqn(protein_data, transform = T, kitwise = F, long_format = F)

  expect_equal(class(normalized_data), class(protein_data))
  expect_equal(dim(normalized_data), dim(protein_data))
  expect_equal(colnames(normalized_data), colnames(protein_data))
  expect_equal(class(apply_protpqn(tibble::as_tibble(protein_data), transform = T,
                                   kitwise = F, long_format = F)),
               class(tibble::as_tibble(protein_data)))

  # With kit info
  kit_info <- data.frame("protein" = proteins,
                         "kit" = rep(c("kitA", "kitB"), each = length(proteins)/2))

  normalized_data <- apply_protpqn(protein_data, transform = T, kitwise = T,
                                   long_format = F, kit_key = kit_info)

  expect_equal(class(normalized_data), class(protein_data))
  expect_equal(dim(normalized_data), dim(protein_data))
  expect_equal(colnames(normalized_data), colnames(protein_data))
  expect_equal(class(apply_protpqn(tibble::as_tibble(protein_data), transform = T,
                                   kitwise = F, long_format = F)),
               class(tibble::as_tibble(protein_data)))

  # Long format test data
  protein_data_long <- data.frame("sample_id" = rep(paste0("s", 1:100), 50),
                                  "protein" = rep(paste0("p", 1:50), each = 100),
                                  "value" = log2(rnorm(100*50, 10, 1)),
                                  "kit" = rep(c("kitA", "kitB"), each = 50*100/2),
                                  "other_column" = "a",
                                  "other_column2" = "b")

  normalized_data_long <- apply_protpqn(protein_data_long, transform = TRUE, kitwise = TRUE, long_format = TRUE)

  expect_equal(class(normalized_data_long), class(protein_data_long))
  expect_equal(dim(normalized_data_long), dim(protein_data_long))
  expect_equal(colnames(normalized_data_long), colnames(protein_data_long))
  expect_equal(class(apply_protpqn(tibble::as_tibble(protein_data_long), transform = TRUE,
                                   kitwise = TRUE, long_format = TRUE)),
               class(tibble::as_tibble(protein_data_long)))
})

test_that("correct errors", {
  # Wide format
  protein_data <- replicate(50, rnorm(100, 10, 1)) |> as.data.frame() |> log2()
  proteins <- paste0("p", 1:50)
  colnames(protein_data) <- proteins
  protein_data$sample_id <- paste0("s", 1:100)

  # No kit info supplied for wide format data when kitwise normalization is specified
  expect_error(apply_protpqn(protein_data, transform = T, kitwise = T, long_format = F),
               "Wide format kitwise normalization requires kit_key")

  kit_info <- data.frame("protein" = proteins,
                         "kit" = rep(c("kitA", "kitB"), each = length(proteins)/2))

  # No sample_id column
  expect_error(apply_protpqn(protein_data |> dplyr::select(-sample_id), transform = T, kitwise = T,
                             long_format = F, kit_key = kit_info),
               "Wide data must contain the column 'sample_id'")
  # Kit info is not a data frame
  expect_error(apply_protpqn(protein_data, transform = T, kitwise = T,
                             long_format = F, kit_key = as.matrix(kit_info)),
               "Kit info .* data frame")
  # Kit info does not contain all proteins
  expect_error(apply_protpqn(protein_data, transform = T, kitwise = T,
                             long_format = F, kit_key = kit_info[1:42, ]),
               "kit_key must contain entries for all")

  # Long format
  protein_data_long <- data.frame("sample_id" = rep(paste0("s", 1:100), 50),
                                  "protein" = rep(paste0("p", 1:50), each = 100),
                                  "value" = log2(rnorm(100*50, 10, 1)),
                                  "kit" = rep(c("kitA", "kitB"), each = 50*100/2),
                                  "other_column" = "a",
                                  "other_column2" = "b")
  # Column names different from requirement or columns not present
  expect_error(apply_protpqn(protein_data_long %>% dplyr::rename(Value = value, Protein = protein),
                             transform = TRUE, kitwise = TRUE, long_format = TRUE),
               "Long data must contain the columns")
  expect_error(apply_protpqn(protein_data_long %>% dplyr::select(-sample_id, -kit), transform = TRUE,
                             kitwise = TRUE, long_format = TRUE),
               "Long data must contain the columns")

  # Protein with multiple measurements per sample (non-unique protein IDs when kitwise is FALSE)
  prot_long <- protein_data_long %>%
    dplyr::filter(protein != "p50")
  prot_long_extra <- prot_long %>%
    dplyr::filter(protein == "p1") %>%
    dplyr::mutate(value = log2(rnorm(nrow(.), 10, 1)),
                  kit = "kitB")
  prot_long <- rbind(prot_long, prot_long_extra)

  expect_error(apply_protpqn(prot_long, transform = TRUE, kitwise = FALSE, long_format = TRUE),
               "unique protein names")
})
