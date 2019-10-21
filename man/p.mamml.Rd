\name{p.mamml}
\alias{p.mamml}
\title{
Compute a combined \emph{p}-value via the model-averaged mean maximized likelihood ratio
}
\description{
The model averaging by mean maximum likelihood (MAMML) test combines likelihood ratio tests and corrects for multiple testing while controlling the weak-sense family-wise error rate in a way that is more powerful than common alternatives including Bonferroni and Simes procedures and possesses theoretical robustness to positive correlations between tests and unequal weights. It is a frequentist analogue to Bayesian model averaging, is theoretically grounded in generalized central limit theorem, and motivates the simpler and better-calibrated harmonic mean \emph{p}-value (HMP) test. The model-averaged mean maximized likelihood (MAMML) is defined as the (possibly weighted) arithmetic mean of the maximized likelihood ratios from a series of likelihood ratio tests comparing mutually exclusive alternative hypotheses with the same nested null hypothesis based on the exact same data.
}
\usage{
p.mamml(R, nu, w = NULL, L = NULL)
}
\arguments{
  \item{R}{
A numeric vector of one or more maximized likelihood ratios. Missing values (NAs) will cause a missing value to be returned.
}
  \item{nu}{
A numeric scalar or vector for the degrees of freedom corresponding to all or each of the maximized likelihood ratios respectively.
}
  \item{w}{
An optional numeric vector of weights that can be interpreted as prior model probabilities for each of the alternative hypotheses represented by the individual \emph{p}-values. The sum of the weights cannot exceed one but may be less than one, which is interpreted as meaning that some \emph{p}-values have been excluded.
}
  \item{L}{
The number of constituent maximized likelihood ratios. If ignored, it defaults to the length of argument \code{R}, with a warning.
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
mamml.stat, hmp.stat, p.hmp
}
\examples{
# For detailed examples type vignette("harmonicmeanp")
nu = 3
R = exp(0.5*rchisq(1000,nu))
mamml.stat(R)
p.mamml(R,nu,L=1000)
}
\keyword{ ~MAMML }
