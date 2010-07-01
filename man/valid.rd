\name{valid}
\alias{valid}

\title{
    Ensure probabilities are valid
}
\description{
  Truncate calculated probabilities into the range [0,1]
}
\usage{
valid(x)
}

\arguments{
  \item{x}{probability}
}
\details{
  Due to floating-point arithmetic, a number that should represent a
    probability can be calculated as being less than zero or greater
    than one.  This function returns a value that is a valid probability.
}
\author{
  Kevin Wright
}
\keyword{ models }
