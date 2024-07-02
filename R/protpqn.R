#' ProtPQN
#'
#' Perform Protein-specific Probabilistic Quotient Normalization on a matrix.
#'
#' @param X Data frame in wide format containing only protein measurements to be normalized. SampleIDs could e.g. be stored as rownames.
#'
#' @return Output is a data frame or matrix in the same format as input with normalized protein levels.
#' @export
#'
#' @examples
#' # Simple example
#' protein_data <- replicate(50, rnorm(100, 10, 1))
#' normalized_data <- protpqn(protein_data)
protpqn <- function(X) {

  stopifnot(all(is.finite(X[!is.na(X)]))) # Ensures finite values
  ref <- apply(X, 2, median, na.rm = T)   # Calculates median for each column
  quo <- t( t(X)/ref )                    # Divides each element by its column-median
  pq <- apply(quo, 1, median, na.rm = T)  # Calculates median for each row from "quo"
  pq[pq == 0] <- NA_real_                 # 0 is replaced by NAs

  # Defines how to handle missing values for correlation calculations
  use <- apply(X, 2, function(x) which(is.na(x))) %>%
    sapply(., function(y) identical(.[[1]], y)) %>% {
      if(all(.)) "complete.obs" else "pairwise.complete.obs"
    }

  # Calculates correlation between row-medians and each column in input
  cor_pq <- cbind(pq, X) %>%
    cor(method= "spearman", use= use) %>%
    .[-1, 1]

  # Antibody specific PQ
  as_pq <- ((pq - 1) %*% t(cor_pq)) + 1   # Adjusts pq-values using correlations values
  out <- X / as_pq                        # Normalize input using adjusted pq-values

  return(out)
}
