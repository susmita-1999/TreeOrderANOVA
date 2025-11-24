\name{TreeMax}
\alias{TreeMax}
\title{Maximum Test Statistic under Tree Order}
\usage{
TreeMax(sample_data, significance_level)
}
\arguments{
  \item{sample_data}{A list-of-lists with treatment-block data.}
  \item{significance_level}{Chosen significance level.}
}
\description{
Computes the MAX test statistic for tree-ordered treatment means in a two-way ANOVA model using bootstrap critical values.
}
\value{
A list containing:
\item{critical_value}{Bootstrap critical value}
\item{test_statistic}{Observed MAX statistic}
\item{decision}{Reject or not reject the null hypothesis}
}
\examples{
## Not run:
# TreeMax(sample_data, 0.05)
## End(Not run)
}
