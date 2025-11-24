\name{TreeMin}
\alias{TreeMin}
\title{Minimum Test Statistic under Tree Order}
\usage{
TreeMin(sample_data, significance_level)
}
\arguments{
  \item{sample_data}{A list-of-lists with data for each treatment and block.}
  \item{significance_level}{Significance level (e.g. 0.05).}
}
\description{
Computes the MIN test statistic for comparing a control treatment to tree-ordered alternatives in a two-way ANOVA layout. Critical values are obtained by bootstrap sampling.
}
\value{
A list containing:
\item{critical_value}{Bootstrap critical value}
\item{test_statistic}{Observed MIN statistic}
\item{decision}{Reject or do not reject the null hypothesis}
}
\examples{
## Not run:
# TreeMin(sample_data, 0.05)
## End(Not run)
}
