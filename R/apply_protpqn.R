#' Apply ProtPQN
#'
#' Apply ProtPQN to wide or long format data with options for data transformation from log and protein kit/panel-wise normalization.
#'
#' @param data_in Data frame containing protein measurement values to be normalized.
#' If the data is in long format, it requires columns named 'sample_id' containing sample identifiers,
#' 'protein' with protein names, 'value' with the values to normalize, and 'kit' with identifiers for kits/panel the proteins are from.
#' If the data is in wide format, it requires a column named 'sample_id' for sample identifiers,
#' and one column for each protein measured containing the values.
#' @param transform Logical. If TRUE, the incoming log-scale data is transformed before normalization and then converted back to log-scale after.
#' @param log_base Numeric of length 1. The base to use when transforming from and to log-scale if transform is TRUE. Default is 2.
#' @param kitwise Logical. If TRUE, normalization is performed kit-wise.
#' @param long_format Logical. If TRUE, the input data is in long format.
#' @param kit_key Data frame containing the kit/panel information for each protein. Required if long_format is FALSE, i.e. data is in wide format.
#'
#' @return A data frame containing normalized protein values in the same format as input.
#' @export
#'
#' @examples
#' # Wide format input example:
#' # Without kit information
#' protein_data <- replicate(50, rnorm(100, 10, 1)) |> as.data.frame() |> log2()
#' proteins <- paste0("p", 1:50)
#' colnames(protein_data) <- proteins
#' protein_data$sample_id <- paste0("s", 1:100)
#' protein_data <- protein_data[, c(ncol(protein_data), 1:(ncol(protein_data)-1))]
#' normalized_data <- apply_protpqn(protein_data,
#'                                  transform = TRUE,
#'                                  kitwise = FALSE,
#'                                  long_format = FALSE)
#' # With kit information
#' kit_info <- data.frame("protein" = proteins,
#'                        "kit" = rep(c("kitA", "kitB"), each = length(proteins)/2))
#'
#' normalized_data <- apply_protpqn(protein_data, transform = TRUE, kitwise = TRUE,
#'                                  long_format = FALSE, kit_key = kit_info)
#'
#' # For long format input:
#' protein_data_long <- data.frame("sample_id" = rep(paste0("s", 1:100), 50),
#'                                 "protein" = rep(paste0("p", 1:50), each = 100),
#'                                 "value" = log2(rnorm(100*50, 10, 1)),
#'                                 "kit" = rep(c("kitA", "kitB"), each = 50*100/2))
#' normalized_data_long <- apply_protpqn(protein_data_long,
#'                                       transform = TRUE,
#'                                       kitwise = TRUE,
#'                                       long_format = TRUE)
apply_protpqn <- function(data_in,
                          transform = TRUE, log_base = 2,
                          kitwise = TRUE, long_format = TRUE, kit_key = NULL) {

  if (long_format && !all(c("value", "protein", "sample_id", "kit") %in% colnames(data_in))) {
    stop("Long data must contain the columns 'sample_id', 'protein', 'value', and 'kit'.")
  }

  # Check that there are no duplicate protein measurements in long-format data (non-unique protein names)
  # May happen with multi-kit data where a protein is included in multiple kits but kitwise is FALSE
  # Wide format already guarantees unique protein names due to uniqueness of column names
  if (long_format & !kitwise) {
    duplicate_proteins <- data_in %>%
      dplyr::count(sample_id, protein) %>%
      dplyr::filter(n > 1L)

    if (nrow(duplicate_proteins) > 0L) {
      stop("Long data without kitwise normalization requires unique protein names (i.e. at most one value per sample_id and protein combination).")
    }
  }

  if (!long_format && !("sample_id" %in% colnames(data_in))) {
    stop("Wide data must contain the column 'sample_id'.")
  }

  if (!long_format && kitwise) {
    if (is.null(kit_key)) {stop("Wide format kitwise normalization requires kit_key.")}
    if (!is.data.frame(kit_key)) {stop("Kit info for wide format data must be given as a data frame.")}
    if (!all(colnames(data_in)[-which(colnames(data_in) == "sample_id")] %in% kit_key$protein)) {stop("kit_key must contain entries for all proteins present in input data.")}
  }

  # Wide input converted to long_format. Keep old column names and class to give output in same format
  if (!long_format) {
    data_cols <- colnames(data_in)
    data_class <- class(data_in)

    data_in <- data_in |>
      tidyr::pivot_longer(cols = -sample_id, names_to = "protein", values_to = "value")

    # Add kit information
    if (!is.null(kit_key)) {
      data_in <- data_in |>
        dplyr::left_join(kit_key, by = "protein")
    } else {
      # If no kit information was provided, add dummy column for kit (not used for normalization as kitwise requires kit input)
      data_in$kit <- "dummy"
    }
  }

  if(kitwise) {
    unique_kits <- unique(data_in$kit) # Get unique values in the "kit" column
    kitwise_data <- vector("list", length(unique_kits)) # Initializes empty list for storage of results

    for (i in seq_along(unique_kits)) {
      kit_name <- unique_kits[i]

      # Filter data for the current kit
      data_kit <- data_in %>%
        dplyr::filter(kit == kit_name) |>
        tidyr::pivot_wider(id_cols = "sample_id", names_from = "protein", values_from = "value") |>
        tibble::column_to_rownames(var = "sample_id")

      data_result <- data_kit |>
        (\(.) { if (transform) log_base^(.) else . })() |> # Transforms the data
        as.matrix() |>
        protpqn() |>
        (\(.) { if (transform) logb(., base = log_base) else . })() |> # Reapplies log
        tibble::as_tibble(rownames = NA) |>
        dplyr::mutate(sample_id = rownames(data_kit)) |>
        tidyr::pivot_longer(cols = -sample_id, names_to = "protein", values_to = "value") |>  # Returns to long_format
        dplyr::mutate(kit = kit_name)
      kitwise_data[[i]] <- data_result # Stores the results

    }

    # Here the kitwise data frames are recombined for the final output
    norm_data <- data.frame(sample_id = character(0), protein = character(0), value = numeric(0), kit = character(0))
    for (df in kitwise_data) {
      norm_data <- dplyr::bind_rows(norm_data, df)
    }
  }

  # Without kitwise normalization
  else {
    key <- data_in |> # Used to add the kit info back to the normalized data frame
      dplyr::select(protein, kit) |>
      dplyr::distinct()

    prep_data <- data_in |>
      tidyr::pivot_wider(id_cols = "sample_id", names_from = "protein", values_from = "value") |>
      tibble::column_to_rownames(var = "sample_id") |>
      (\(.) { if (transform) log_base^(.) else . })() |>
      as.matrix()

    norm_data <- prep_data |>
      protpqn() |>
      (\(.) { if (transform) logb(., base = log_base) else . })() |>
      tibble::as_tibble(rownames = NA) |>
      dplyr::mutate(sample_id = rownames(prep_data)) |>
      tidyr::pivot_longer(cols = -sample_id, names_to = "protein", values_to = "value") |>
      dplyr::left_join(key, by = "protein")
  }

  # Give out put in same format and class as input
  if (long_format) {
    value_pos <- which(colnames(data_in) == "value")
    norm_data <- data_in %>%
      dplyr::select(-value) %>%
      dplyr::left_join(norm_data, by = c("sample_id", "protein", "kit")) %>%
      dplyr::relocate(value, .after = value_pos - 1)
  } else {
    norm_data <- norm_data |>
      tidyr::pivot_wider(id_cols = "sample_id", names_from = "protein", values_from = "value")
    norm_data <- norm_data[, data_cols]

    if (!"tbl" %in% data_class) {
      norm_data <- as.data.frame(norm_data)
    }
  }

  return(norm_data)
}
