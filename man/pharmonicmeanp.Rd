\name{Harmonic Mean p-Value Asymptotic Distribution}
\alias{dharmonicmeanp}
\alias{pharmonicmeanp}
\alias{qharmonicmeanp}
\alias{rharmonicmeanp}
\title{
The Harmonic Mean p-Value Asymptotic Distribution
}
\description{
Density, distribution function, quantile function and random number generation for the harmonic mean of \emph{L} \emph{p}-values under their null hypotheses, i.e. the harmonic mean of \emph{L} standard uniform random variables, assuming \emph{L} is large.
}
\usage{
dharmonicmeanp(x, L, log=FALSE)
pharmonicmeanp(x, L, log=FALSE, lower.tail=TRUE)
qharmonicmeanp(p, L, log=FALSE, lower.tail=TRUE)
rharmonicmeanp(n, L)
}
\arguments{
  \item{x}{
The value or vector of values of the harmonic mean \emph{p}-value, for example calculated from data using function \code{hmp.stat}.
}
  \item{L}{
The number of constituent \emph{p}-values used in calculating each value of x.
}
  \item{log}{
If true the log probability is output.
}
  \item{lower.tail}{
If true (the default) the lower tail probability is returned. Otherwise the upper tail probability.
}
  \item{p}{
The value or vector of values, between 0 and 1, of the probability specifying the quantile for which to return the harmonic mean \emph{p}-value.
}
  \item{n}{
The number of values to simulate.
}
}
\value{
\code{dharmonicmeanp} produces the density, \code{pharmonicmeanp} the tail probability, \code{qharmonicmeanp} the quantile and \code{rharmonicmeanp} random variates for the harmonic mean of \emph{L} \emph{p}-values when their null hypotheses are true.

Use \code{qharmonicmeanp(alpha,L)} to calculate \eqn{\alpha_L}, the `harmonic mean \emph{p}-value threshold', as in Table 1 of Wilson (2019, corrected), where \code{L} is the total number of \emph{p}-values under consideration and \code{alpha} is the intended strong-sense familywise error rate.

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
# For a detailed tutorial type vignette("harmonicmeanp")
# Example: simulate from a non-uniform distribution mildly enriched for small \emph{p}-values. 
# Compare the significance of the combined p-value for Bonferroni, Benjamini-Hochberg (i.e. Simes), 
# HMP and (equivalently) MAMML with 2 degrees of freedom.
L = 1000
p = rbeta(L,1/1.5,1)
min(p.adjust(p,"bonferroni"))
min(p.adjust(p,"BH"))
x = hmp.stat(p)
pharmonicmeanp(x,length(p))
p.hmp(p,L=L)
p.mamml(1/p,2,L=L)

# Compute critical values for the HMP from asymptotic theory and compare to direct simulations
L = 100
alpha = 0.05
(hmp.crit = qharmonicmeanp(alpha,L))
nsim = 100000
p.direct = matrix(runif(L*nsim),nsim,L)
hmp.direct = apply(p.direct,1,hmp.stat)
(hmp.crit.sim = quantile(hmp.direct,alpha))

# Compare HMP of p-values simulated directly, and via the asymptotic distribution, 
# to the asymptotic density
L = 30
nsim = 10000
p.direct = matrix(runif(L*nsim),nsim,L)
hmp.direct = apply(p.direct,1,hmp.stat)
hmp.asympt = rharmonicmeanp(nsim,L)
h = hist(hmp.direct,60,col="green3",prob=TRUE,main="Distributions of harmonic mean p-values")
hist(hmp.asympt,c(-Inf,h$breaks,Inf),col="yellow2",prob=TRUE,add=TRUE)
hist(hmp.direct,60,col=NULL,prob=TRUE,add=TRUE)
curve(dharmonicmeanp(x,L),lwd=2,col="red3",add=TRUE)
legend("topright",c("Direct simulation","Asymptotic simulation","Asymptotic density"),
  fill=c("green3","yellow2",NA),col=c(NA,NA,"red3"),lwd=c(NA,NA,2),bty="n",border=c(1,1,NA))
}
\keyword{ ~hmp }
