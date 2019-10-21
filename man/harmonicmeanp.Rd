\name{harmonicmeanp}
\alias{harmonicmeanp}
\alias{hmp}
\docType{package}
\title{
Harmonic mean p-values and model averaging by mean maximum likelihood
}
\description{
The harmonic mean \emph{p}-value (HMP) test combines \emph{p}-values and corrects for multiple testing while controlling the strong-sense family-wise error rate. It is more powerful than common alternatives including Bonferroni and Simes procedures when combining large proportions of all the \emph{p}-values, at the cost of slightly lower power when combining small proportions of all the \emph{p}-values. It is more stringent than controlling the false discovery rate, and possesses theoretical robustness to positive correlations between tests and unequal weights. It is a multi-level test in the sense that a superset of one or more significant tests is certain to be significant and conversely when the superset is non-significant, the constituent tests are certain to be non-significant. It is based on MAMML (model averaging by mean maximum likelihood), a frequentist analogue to Bayesian model averaging, and is theoretically grounded in generalized central limit theorem.
}
\details{
\tabular{ll}{
Package: \tab harmonicmeanp\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2019-08-17\cr
License: \tab GPL-3\cr
}
The key function is p.hmp for combining p-values using the HMP. Type vignette("harmonicmeanp") for detailed examples.
}
\author{
Daniel J. Wilson

Maintainer: Daniel Wilson <hmp.R.package@gmail.com>
}
\references{
Daniel J. Wilson (2019) The harmonic mean \emph{p}-value for combining dependent tests. \emph{Proceedings of the National Academy of Sciences USA} 116: 1195-1200.
}
\seealso{
Package FMStable
}
\examples{
# For detailed examples type vignette("harmonicmeanp")
# Example: simulate from a non-uniform distribution mildly enriched for small \emph{p}-values. 
# Compare the significance of the combined p-value for Bonferroni, Benjamini-Hochberg (i.e. Simes), 
# HMP and (equivalently) MAMML with 2 degrees of freedom.
p = rbeta(1000,1/1.5,1)
min(p.adjust(p,"bonferroni"))
min(p.adjust(p,"BH"))
p.hmp(p,L=1000)
p.mamml(1/p,2,L=1000)
}
