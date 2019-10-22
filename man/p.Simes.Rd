\name{Multilevel Simes test}
\alias{p.Simes}
\title{
Multilevel Simes test
}
\description{
Compute a combined \emph{p}-value via the multilevel Simes test. This function is an alternative to \code{p.hmp} which calculates the asymptotically exact harmonic mean \emph{p}-value (HMP). Simes' test is more robust to certain types of dependence than the HMP, and might be preferred in some cases.}
\usage{
p.Simes(p, w = NULL, L = NULL, w.sum.tolerance = 1e-6, multilevel = TRUE)
}
\arguments{
  \item{p}{
A numeric vector of one or more \emph{p}-values to be combined. Missing values (NAs) will cause a missing value to be returned.
}
  \item{w}{
An optional numeric vector of weights that can be interpreted as incorporating information about prior model probabilities and power of the tests represented by the individual \emph{p}-values. The sum of the weights cannot exceed one but may be less than one, which is interpreted as meaning that some of the \code{L} \emph{p}-values have been excluded.
}
  \item{L}{
The number of constituent \emph{p}-values. This is mandatory when using \code{p.Simes} as part of a multilevel test procedure and needs to be set equal to the total number of individual \emph{p}-values investigated, which may be (much) larger than the length of the argument \code{p}. If ignored, it defaults to the length of argument \code{p}, with a warning.
}
  \item{w.sum.tolerance}{
Tolerance for checking that the weights do not exceed 1.
}
  \item{multilevel}{
Logical, indicating whether the test is part of a multilevel procedure involving a total of \code{L} \emph{p}-values intended to control the strong-sense familywise error rate (\code{TRUE}, default), or whether a stand-alone test is to be produced (\code{FALSE}), in which case \code{L} is ignored.
}
}
\value{
A \emph{p}-value based on Simes' procedure is returned. This \emph{p}-value can be compared against threshold \code{sum(w)*alpha} to control the strong-sense familywise error rate at level \code{alpha}. A test based on this asymptotically exact harmonic mean \emph{p}-value is equivalent to comparison of the `raw' harmonic mean \emph{p}-value calculated by \code{hmp.stat(p,w)} to a `harmonic mean \emph{p}-value threshold' \code{sum(w)*qharmonicmeanp(alpha,L)}.

Unlike the Benjamini-Hochberg (BH) procedure, which uses Simes' procedure to control the false discovery rate, the multilevel Simes test controls the family-wise error rate. Whereas the BH procedure tests the significance of \emph{individual} hypotheses, the multilevel Simes procedure is a combined test of whether a \emph{group} of hypotheses are statistically significant. If the group comprises a single \emph{p}-value, the procedure defaults to Bonferroni correction. Like the HMP, the multilevel Simes procedure can be used to find the smallest groups of hypotheses that are significant.
}
\references{
R. J. Simes (1986) An improved Bonferroni procedure for multiple tests of significance. \emph{Biometrika} 73:751–754.
D. J. Wilson (2019) The harmonic mean \emph{p}-value for combining dependent tests. \emph{Proceedings of the National Academy of Sciences USA} 116: 1195-1200.
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
p = rbeta(L,1/1.8,1)
min(p.adjust(p,"bonferroni"))
min(p.adjust(p,"BH"))
p.hmp(p,L=L)
p.mamml(1/p,2,L=L)
p.Simes(p,L=L)

# Multilevel test: find significant subsets of the 1000 p-values by comparing overlapping
# subsets of size 100 and 10 in addition to the top-level test of all 1000, while maintaining
# the strong-sense family-wise error rate.
p.100 = sapply(seq(1,L-100,by=10),function(beg) p.hmp(p[beg+0:99],L=L))
s.100 = sapply(seq(1,L-100,by=10),function(beg) p.Simes(p[beg+0:99],L=L))
plot(-log10(p.100),xlab="Test index",main="Moving average",
  ylab="Asymptotically exact combined -log10 p-value",type="o")
lines(-log10(s.100),type="o",col="yellow3")
# The appropriate threshold is alpha (e.g. 0.05) times the proportion of tests in each p.100
abline(h=-log10(0.05*100/1000),col=2,lty=2)

p.10 = sapply(seq(1,L-10,by=5),function(beg) p.hmp(p[beg+0:9],L=L))
s.10 = sapply(seq(1,L-10,by=5),function(beg) p.Simes(p[beg+0:9],L=L))
plot(-log10(p.10),xlab="Test index",main="Moving average",
  ylab="Asymptotically exact combined -log10 p-value",type="o")
lines(-log10(s.10),type="o",col="yellow3")
# The appropriate threshold is alpha (e.g. 0.05) times the proportion of tests in each p.10
abline(h=-log10(0.05*10/1000),col=2,lty=2)
}
\keyword{ ~Simes }
