
# ProtPQN

<!-- badges: start -->
<!-- badges: end -->

Apply Protein-specific Probabilistic Quotient Normalization (ProtPQN) to affinity proteomics data from liquid biopsies. Use on any wide-format matrix or data frame containing protein measurements, or use the wrapper function to apply on wide or long data with options for transforming back from log-scale data before normalization and normalization in a kit/panel-wise manner. 

The method is based on Probabilistic Quotient Normalization (Dieterle et al., 2006) and aims to reduce sample-to-sample variation introduced by uneven handling of samples, while adjusting the normalization factor for each protein to take into account differences in the affinities of the affinity reagents. ProtPQN has been described by [Dodig-Crnković et al. (2020)](https://doi.org/10.1016%2Fj.ebiom.2020.102854) and [Fredolini et al. (2024)](https://doi.org/10.1038/s43856-024-00480-4), and implemented for use with Suspension Bead Array data in our [BAf package](https://github.com/Schwenk-Lab/BAf-R_package) (under the name 'AbsPQN'). This package simplifies the application of the method for any affinity proteomics data generated from liquid biopsies.

## Installation

You can install the development version of ProtPQN like so:

``` r
# install.packages("devtools")
devtools::install_github("Schwenk-Lab/ProtPQN")
```

## Example

The following example shows a simple case with wide data in a matrix, and another case with more complicated data that has been log2-transformed and contains multiple protein kits.

``` r
library(ProtPQN)
# Simple example on wide matrix
protein_data <- replicate(50, rnorm(100, 10, 1))
normalized_data <- protpqn(protein_data)

# Data transformation and kit-wise normalization, here on long data
protein_data_long <- data.frame("sample_id" = rep(paste0("s", 1:100), 50),
                                "protein" = rep(paste0("p", 1:50), each = 100),
                                "value" = log2(rnorm(100*50, 10, 1)),
                                "kit" = rep(c("kitA", "kitB"), each = 50*100/2))
normalized_data_long <- apply_protpqn(protein_data_long,
                                      transform = TRUE,
                                      kitwise = TRUE,
                                      long_format = TRUE)
```

