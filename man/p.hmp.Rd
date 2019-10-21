\name{Asymptotically Exact Harmonic Mean p-Value}
\alias{p.hmp}
\title{
Asymptotically Exact Harmonic Mean p-Value
}
\description{
Compute a combined \emph{p}-value via the asymptotically exact harmonic mean \emph{p}-value. The harmonic mean \emph{p}-value (HMP) test combines \emph{p}-values and corrects for multiple testing while controlling the strong-sense family-wise error rate. It is more powerful than common alternatives including Bonferroni and Simes procedures when combining large proportions of all the \emph{p}-values, at the cost of slightly lower power when combining small proportions of all the \emph{p}-values. It is more stringent than controlling the false discovery rate, and possesses theoretical robustness to positive correlations between tests and unequal weights. It is a multi-level test in the sense that a superset of one or more significant tests is certain to be significant and conversely when the superset is non-significant, the constituent tests are certain to be non-significant. It is based on MAMML (model averaging by mean maximum likelihood), a frequentist analogue to Bayesian model averaging, and is theoretically grounded in generalized central limit theorem.
}
\usage{
p.hmp(p, w = NULL, L = NULL, w.sum.tolerance = 1e-6, multilevel = TRUE)
}
\arguments{
  \item{p}{
A numeric vector of one or more \emph{p}-values to be combined. Missing values (NAs) will cause a missing value to be returned.
}
  \item{w}{
An optional numeric vector of weights that can be interpreted as incorporating information about prior model probabilities and power of the tests represented by the individual \emph{p}-values. The sum of the weights cannot exceed one but may be less than one, which is interpreted as meaning that some of the \code{L} \emph{p}-values have been excluded.
}
  \item{L}{
The number of constituent \emph{p}-values. This is mandatory when using \code{p.hmp} as part of a multilevel test procedure and needs to be set equal to the total number of individual \emph{p}-values investigated, which may be (much) larger than the length of the argument \code{p}. If ignored, it defaults to the length of argument \code{p}, with a warning.
}
  \item{w.sum.tolerance}{
Tolerance for checking that the weights do not exceed 1.
}
  \item{multilevel}{
Logical, indicating whether the test is part of a multilevel procedure involving a total of \code{L} \emph{p}-values intended to control the strong-sense familywise error rate (\code{TRUE}, default), or whether a stand-alone test is to be produced (\code{FALSE}), in which case \code{L} is ignored.
}
}
\value{
When \code{multilevel==TRUE} an asymptotically exact combined \emph{p}-value is returned, equivalent to \code{sum(w)*pharmonicmeanp(hmp.stat(p,w)/sum(w),L)}. This \emph{p}-value can be compared against threshold \code{sum(w)*alpha} to control the strong-sense familywise error rate at level \code{alpha}. A test based on this asymptotically exact harmonic mean \emph{p}-value is equivalent to comparison of the `raw' harmonic mean \emph{p}-value calculated by \code{hmp.stat(p,w)} to a `harmonic mean \emph{p}-value threshold' \code{sum(w)*qharmonicmeanp(alpha,L)}.

When \code{multilevel==FALSE} an asymptotically exact combined \emph{p}-value is returned, equivalent to \code{pharmonicmeanp(hmp.stat(p,w),length(p))}. \code{L} is ignored and \code{w} is normalized to sum to 1. This \emph{p}-value can be compared against threshold \code{alpha} to control the per-test error rate at level \code{alpha}. A test based on this asymptotically exact harmonic mean \emph{p}-value is equivalent to comparison of the `raw' harmonic mean \emph{p}-value calculated by \code{hmp.stat(p,w)} to a `harmonic mean \emph{p}-value threshold' \code{qharmonicmeanp(alpha,L)}.

}
\references{
Daniel J. Wilson (2019) The harmonic mean \emph{p}-value for combining dependent tests. \emph{Proceedings of the National Academy of Sciences USA} 116: 1195-1200.
}
\author{
Daniel J. Wilson
}

\seealso{
hmp.stat
}
\examples{
# For detailed examples type vignette("harmonicmeanp")
# Example: simulate from a non-uniform distribution mildly enriched for small \emph{p}-values. 
# Compare the significance of the combined p-value for Bonferroni, Benjamini-Hochberg (i.e. Simes), 
# HMP and (equivalently) MAMML with 2 degrees of freedom.
L = 1000
p = rbeta(L,1/1.8,1)
min(p.adjust(p,"bonferroni"))
min(p.adjust(p,"BH"))
p.hmp(p,L=L)
p.mamml(1/p,2,L=L)

# Multilevel test: find significant subsets of the 1000 p-values by comparing overlapping
# subsets of size 100 and 10 in addition to the top-level test of all 1000, while maintaining
# the strong-sense family-wise error rate.
p.100 = sapply(seq(1,L-100,by=10),function(beg) p.hmp(p[beg+0:99],L=L))
plot(-log10(p.100),xlab="Test index",main="Moving average",
  ylab="Asymptotically exact combined -log10 p-value",type="o")
# The appropriate threshold is alpha (e.g. 0.05) times the proportion of tests in each p.100
abline(h=-log10(0.05*100/1000),col=2,lty=2)

p.10 = sapply(seq(1,L-10,by=5),function(beg) p.hmp(p[beg+0:9],L=L))
plot(-log10(p.10),xlab="Test index",main="Moving average",
  ylab="Asymptotically exact combined -log10 p-value",type="o")
# The appropriate threshold is alpha (e.g. 0.05) times the proportion of tests in each p.10
abline(h=-log10(0.05*10/1000),col=2,lty=2)
}
\keyword{ ~hmp }
