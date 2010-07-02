\name{stack2}
\alias{stack2}
\alias{stack3}
\alias{stack2Excel}
\alias{stack3Excel}

\title{Estimating presence of stacked genes}

\description{
  Assuming qualitative tests are performed on n pools
  of m seeds, estimate the percent of seeds with single 
  genetic traits and the percentage of seeds with stacked
  genetic traits.
}

\usage{
stack2(n, m, nA, nB, nAB, existAB="Yes", fpr=0, fnr=0, check=TRUE)
stack3(n, m, nA, nB, nC, nAB, nAC, nBC, nABC,
       existAB="Yes", existAC="Yes", existBC="Yes", existABC="Yes",
       fpr=0, fnr=0, check=TRUE)
stack2Excel(...)
stack3Excel(...)

}

\arguments{
  \item{n}{the number of pools}
  \item{m}{the number of seeds in each pool}
  \item{nA}{the number of positive pools for event A only}
  \item{nB}{the number of positive pools for event B only}
  \item{nAB}{the number of positive pools for both A and B}
  \item{nC}{the number of positive pools for event C only}  
  \item{nAC}{the number of positive pools for both A and C}
  \item{nBC}{the number of positive pools for both B and C}
  \item{nABC}{the number of positive pools for both A and B and C}
  \item{existAB}{do seeds with a stacked event 'AB' exist?}
  \item{existAC}{do seeds with a stacked event 'AC' exist?}
  \item{existBC}{do seeds with a stacked event 'BC' exist?}
  \item{existABC}{do seeds with a stacked event 'ABC' exist?}
  \item{fpr}{false positive rate (proportion) for detecting GM events}
  \item{fnr}{false negative rate (proportion) for detecting GM events}
  \item{check}{Should simple checks be performed?  Defaults to TRUE}
  \item{...}{Other arguments passed}
}

\value{
  A data frame with the estimated proportion of seeds for each event,
  the observed and expected number of positive pools,
  and whether or not each event can exist.
}

\details{
  The 'stack2Excel' and 'stack3Excel' functions are simple wrappers that
  are intended to be called from Excel and should not issue any warnings.
}

\author{Kevin Wright, Jean-Louis Laffont}

\examples{
stack2(10, 300, 0, 1, 2)
stack3(20,150, 2,2,2,2,2,2,3, existAB="no", fnr=.02, fpr=.02)
}

\keyword{models}
