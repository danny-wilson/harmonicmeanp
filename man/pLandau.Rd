\name{Landau Distribution}
\alias{dLandau}
\alias{pLandau}
\alias{qLandau}
\alias{rLandau}
\title{
The Landau Distribution
}
\description{
Density, distribution function, quantile function and random number generation for the Landau distribution with location parameter mu and scale parameter sigma.
}
\usage{
dLandau(x, mu=log(pi/2), sigma=pi/2, log=FALSE)
pLandau(x, mu=log(pi/2), sigma=pi/2, log=FALSE, lower.tail=TRUE)
qLandau(p, mu=log(pi/2), sigma=pi/2, log=FALSE, lower.tail=TRUE)
rLandau(n, mu=log(pi/2), sigma=pi/2)
}
\arguments{
  \item{x}{
The value or vector of values of the Landau-distributed random variable.
}
  \item{mu}{
The location parameter of the Landau distribution. Defaults to log(pi/2) to give Landau's original distribution.
}
  \item{sigma}{
The scale parameter of the Landau distribution. Defaults to pi/2 to give Landau's original distribution.
}
  \item{log}{
If true the log probability is output.
}
  \item{lower.tail}{
If true (the default) the lower tail probability is returned. Otherwise the upper tail probability.
}
  \item{p}{
The value or vector of values, between 0 and 1, of the probability specifying the quantile for which to return the Landau random variable \code{x}.
}
  \item{n}{
The number of values to simulate.
}
}
\details{
The density of the Landau distribution can be written
\deqn{f(x)=\frac{1}{\pi\,\sigma}\int_0^\infty \exp\left(-t\frac{(x-\mu)}{\sigma}-\frac{2}{\pi}t\log(t)\right)\,\sin\left(2t\right)\textrm{d}t}{f(x)=\frac{1}{\pi\,\sigma}\int_0^\infty \exp\left(-t\frac{(x-\mu)}{\sigma}-\frac{2}{\pi}t\log(t)\right)\,\sin\left(2t\right)\textrm{d}t}
}
\value{
\code{dLandau} produces the density, \code{pLandau} the tail probability, \code{qLandau} the quantile and \code{rLandau} random variates for the Landau distribution.
}
\references{
Landau LD (1944) On the energy loss of fast particles by ionization. \emph{J Phys USSR} 8:201-205.

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
# Example: simulate from a non-uniform distribution mildly enriched for small \emph{p}-values. 
# Compare the significance of the combined p-value for Bonferroni, Benjamini-Hochberg (i.e. Simes), 
# HMP and (equivalently) MAMML with 2 degrees of freedom.
L = 1000
p = rbeta(L,1/1.5,1)
min(p.adjust(p,"bonferroni"))
min(p.adjust(p,"BH"))
x = hmp.stat(p)
pLandau(1/x,log(length(p))+(1 + digamma(1) - log(2/pi)),pi/2,lower.tail=FALSE)
p.hmp(p,L=L)
p.mamml(1/p,2,L=L)
}
\keyword{ ~hmp }
