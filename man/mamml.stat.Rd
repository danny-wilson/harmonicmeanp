\name{Compute the Model-Averaged Mean Maximized Likelihood}
\alias{mamml.stat}
\title{
Compute the Model-Averaged Mean Maximized Likelihood
}
\description{
The model-averaged mean maximized likelihood (MAMML) is defined as the (possibly weighted) arithmetic mean of the maximized likelihood ratios from a series of likelihood ratio tests comparing mutually exclusive alternative hypotheses with the same nested null hypothesis based on the exact same data.
}
\usage{
mamml.stat(R, w = NULL)
}
\arguments{
  \item{R}{
A numeric vector of one or more maximized likelihood ratios. Missing values (NAs) will cause a missing value to be returned.
}
  \item{w}{
An optional numeric vector of weights that can be interpreted as prior model probabilities for each of the alternative hypotheses represented by the individual \emph{p}-values. The sum of the weights cannot exceed one but may be less than one, which is interpreted as meaning that some \emph{p}-values have been excluded.
}
}
\value{
The model-averaged mean maximized likelihood ratio is returned.
}
\references{
Daniel J. Wilson (2019) The harmonic mean \emph{p}-value for combining dependent tests. \emph{Proceedings of the National Academy of Sciences USA} 116: 1195-1200.
}
\author{
Daniel J. Wilson
}

\seealso{
p.mamml
}
\examples{
# For detailed examples type vignette("harmonicmeanp")
nu = 3
R = exp(0.5*rchisq(1000,nu))
mamml.stat(R)
p.mamml(R,nu,L=1000)
}
\keyword{ ~MAMML }
