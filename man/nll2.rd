\name{nll2}
\alias{nll2}
\alias{nll3}
\title{
  Negative loglikelihood
}
\description{
  Calculate the negative loglikelihood for a two-way or three-way
  stacked event.
}
\usage{
nll2(theta, nPos, m, K)
nll3(theta, nPos, m, K)
}

\arguments{
  \item{theta}{parameters to be optimized}
  \item{nPos}{vector of number of positive pools}
  \item{m}{number of seats per pool}
  \item{K}{coefficient matrix for false positive and negative rate
    adjustments}
}

\details{
  This is not intended to be a user-callable function.
}

\author{
  Kevin Wright, Jean-Louis Laffont
}

\keyword{models}
