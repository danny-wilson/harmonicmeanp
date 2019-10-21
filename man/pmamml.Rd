\name{Model-Averaged Mean Maximized Likelihood Ratio Asymptotic Distribution}
\alias{dmamml}
\alias{pmamml}
\alias{qmamml}
\alias{rmamml}
\title{
Model-Averaged Mean Maximized Likelihood Ratio Asymptotic Distribution
}
\description{
Density, distribution function, quantile function and random number generation for the mean of \emph{L} maximized likelihood ratios under their null hypotheses, i.e. the mean of \emph{L} Pareto(1,1) variables, assuming \emph{L} is large.
}
\usage{
dmamml(x, L, df, log=FALSE)
pmamml(x, L, df, log=FALSE, lower.tail=TRUE)
qmamml(p, L, df, log=FALSE, lower.tail=TRUE, xmin=1+1e-12, xmax=1e12)
rmamml(n, L, df)
}
\arguments{
  \item{x}{
The value or vector of values of the mean maximized likelihood ratios, for example calculated from data using function \code{mamml.stat}.
}
  \item{L}{
The number of constituent likelihood ratios used in calculating each value of x.
}
  \item{df}{
Degrees of freedom of the constituent likelihood ratios, assumed all equal for each value of x.
}
  \item{log}{
If true the log probability is output.
}
  \item{lower.tail}{
If true (the default) the lower tail probability is returned. Otherwise the upper tail probability.
}
  \item{p}{
The value or vector of values, between 0 and 1, of the probability specifying the quantile for which to return the mean maximized likelihood ratio.
}
  \item{xmin, xmax}{
The range of values of the mean maximized likelihood ratio over which to search for the quantile specified by \code{p}, used by the numerical root finding algorithm. Must bracket the solution.
}
  \item{n}{
The number of values to simulate.
}
}
\value{
\code{dmamml} produces the density, \code{pmamml} the tail probability, \code{qmamml} the quantile and \code{rmamml} random variates for the model-averaged mean maximized likelihood ratios when their null hypotheses are true.

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
# Example: simulate from a gamma distribution mildly enriched for large likelihood ratios.
# Compare the significance of the combined p-value for Bonferroni, Benjamini-Hochberg (i.e. Simes), 
# MAMML with 4 degrees of freedom.
L = 100
df = 4
enrich = 1.5
y = exp(0.5*rgamma(L,shape=df/2,rate=1/2/enrich))
p = pchisq(2*log(y),df,lower.tail=FALSE)
min(p.adjust(p,"bonferroni"))
min(p.adjust(p,"BH"))
(x = mamml.stat(y))
pmamml(x,L,df,lower.tail=FALSE)
p.mamml(y,df,L=L)
# Compare to the HMP - MAMML may act more conservatively because of adjustments for its
# poorer asymptotic approximation
x = hmp.stat(p)
p.hmp(p,L=L)

# Compute critical values for MAMML from asymptotic theory and compare to direct simulations
L = 100
df = 1
alpha = 0.05
(mamml.crit = qmamml(1-alpha,L,df))
nsim = 10000
y.direct = matrix(exp(0.5*rchisq(L*nsim,df)),nsim,L)
mamml.direct = apply(y.direct,1,mamml.stat)
# mamml.crit may be more conservative than mamml.crit.sim because of adjustments for its
# poorer asymptotic approximation
(mamml.crit.sim = quantile(mamml.direct,1-alpha))

# Compare MAMML simulated directly, and via the asymptotic distribution, to the asymptotic density
# Works best for df = 2
L = 30
df = 3
nsim = 1000
y.direct = matrix(exp(0.5*rchisq(L*nsim,df)),nsim,L)
mamml.direct = apply(y.direct,1,mamml.stat)
xmax = quantile(mamml.direct,.95)
h = hist(mamml.direct,c(-Inf,seq(0,xmax,len=60),Inf),col="green3",prob=TRUE,
  main="Distributions of MAMML",xlim=c(1,xmax))
# Slow because rmamml calls qmamml which uses numerical root finding
# mamml.asympt = rmamml(nsim,L,df)
# hist(mamml.asympt,c(-Inf,seq(0,xmax,len=60),Inf),col="yellow2",prob=TRUE,add=TRUE)
hist(mamml.direct,c(-Inf,seq(0,xmax,len=60),Inf),col=NULL,prob=TRUE,add=TRUE)
curve(dmamml(x,L,df),lwd=2,col="red3",add=TRUE)
legend("topright",c("Direct simulation","Asymptotic density"),
  fill=c("green3",NA),col=c(NA,"red3"),lwd=c(NA,2),bty="n",border=c(1,NA))
}
\keyword{ ~hmp }
