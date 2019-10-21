\name{Compute the Harmonic Mean p-Value}
\alias{hmp.stat}
\title{
Compute the Harmonic Mean \emph{p}-Value
}
\description{
The harmonic mean \emph{p}-value (HMP) is defined as the inverse of the (possibly weighted) arithmetic mean of the inverse \emph{p}-values. When the HMP is small (e.g. less than 0.05), it is approximately well-calibrated, meaning that it can be directly interpreted. However, the function p.hmp calculates an asymptotically exact \emph{p}-value from the HMP and is preferred.
}
\usage{
hmp.stat(p, w = NULL)
}
\arguments{
  \item{p}{
A numeric vector of one or more \emph{p}-values. Missing values (NAs) will cause a missing value to be returned.
}
  \item{w}{
An optional numeric vector of weights that can be interpreted as prior model probabilities for each of the alternative hypotheses represented by the individual \emph{p}-values. The sum of the weights cannot exceed one but may be less than one, which is interpreted as meaning that some \emph{p}-values have been excluded.
}
}
\value{
The harmonic mean \emph{p}-value is returned.
}
\references{
Daniel J. Wilson (2019) The harmonic mean \emph{p}-value for combining dependent tests. \emph{Proceedings of the National Academy of Sciences USA} 116: 1195-1200.
}
\author{
Daniel J. Wilson
}

\seealso{
p.hmp
}
\examples{
# For detailed examples type vignette("harmonicmeanp")
p = rbeta(1000,1/1.5,1)
hmp.stat(p)
p.hmp(p,L=1000)
}
\keyword{ ~hmp }
