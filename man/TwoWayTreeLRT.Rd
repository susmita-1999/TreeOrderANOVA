\name{TwoWayTreeLRT}
\alias{TwoWayTreeLRT}
\title{Two-Way Likelihood Ratio Test Under Tree-Order Restriction}
\description{
Performs the Likelihood Ratio Test (LRT) for testing equality of row
effects in a two-way ANOVA model under a tree-ordered alternative.
}
\usage{
TwoWayTreeLRT(sample_data, significance_level)
}
\arguments{
  \item{sample_data}{A list-of-lists, where
  \code{sample_data[[i]][[j]]} contains observations in cell (i, j).}
  \item{significance_level}{Significance level (default: 0.05).}
}
\value{
A list containing:
  \item{critical_value}{Bootstrap critical value.}
  \item{test_statistic}{Observed LRT statistic.}
  \item{decision}{Reject or not reject the null hypothesis.}
}
\examples{
# sample_data = list(
#   list(rnorm(20), rnorm(20)),
#   list(rnorm(20,1), rnorm(20,1))
# )
# TwoWayTreeLRT(sample_data, 0.05)
}
