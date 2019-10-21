library(FMStable)
hmp.stat = function(p, w = NULL) {
    p = as.numeric(p)
    if (is.null(w)) return(c(hmp.stat = 1/mean(1/p)))
    return(c(hmp.stat = sum(w)/sum(w/p)))
}
p.hmp = function(p, w = NULL, L = NULL, w.sum.tolerance = 1e-6, multilevel = TRUE) {
    if(is.null(L) & multilevel) {
        warning("L not specified: for multilevel testing set L to the total number of individual p-values")
        L = length(p)
    }
    if(length(p) == 0) return(NA)
    if(length(p) > L) warning("The number of p-values cannot exceed L")
    if(is.null(w)) {
        w = rep(1/L,length(p))
    } else if(any(w<0)) {
        stop("No weights can be negative")
    }
    w.sum = sum(w)
    if(w.sum>1+w.sum.tolerance) {
        stop("Weights cannot exceed 1")
    }
    HMP = hmp.stat(p, w)
    O.874 = 1 + digamma(1) - log(2/pi)
    if(multilevel) {
        return(c(p.hmp = w.sum*pEstable(w.sum/HMP, setParam(alpha = 1, location = (log(L) + O.874), logscale = log(pi/2), pm = 0), lower.tail = FALSE)))
    }
    return(c(p.hmp = pEstable(1/HMP, setParam(alpha = 1, location = (log(length(p)) + O.874), logscale = log(pi/2), pm = 0), lower.tail = FALSE)))
}
mamml.stat = function(R, w = NULL) {
    R = as.numeric(R)
    if(any(R<1)) stop("Maximized likelihood ratios cannot be less than one")
    if(is.null(w)) return(c("mamml.stat"=mean(R)))
    w = w/sum(w)
    return(c("mamml.stat"=sum(w*R)))
}
p.mamml = function(R, nu, w = NULL, L = NULL) {
    if(is.null(L)) {
        warning("L not specified, assuming L = length(R)")
        L = length(R)
    }
    if(length(nu)!=1 & length(nu)!=length(R)) stop("Degrees of freedom (nu) must have length one or length of R")
    Rbar = mamml.stat(R, w)
    nu = as.numeric(nu)
    if(any(nu<=0)) stop("Degrees of freedom (nu) must be positive")
    nu.max = max(nu)
    if(nu.max<2) {
        c = pgamma(log(Rbar),nu.max/2,1,lower.tail=FALSE)*Rbar
    } else {
        c = pgamma(log(length(R)*Rbar),nu.max/2,1,lower.tail=FALSE)*length(R)*Rbar
    }
    O.874 = 1+digamma(1)-log(2/pi)
    return(c("p.mamml"=pEstable(Rbar,setParam(alpha=1,location=c*(log(L)+O.874),logscale=log(pi/2*c),pm=0),lower.tail=FALSE)))
}
dLandau = Vectorize(function(x,mu=log(pi/2),sigma=pi/2,log=FALSE) {
    param = setParam(alpha=1, location=mu, logscale=log(sigma), pm=0)
    return(dEstable(x,param,log=log))
})
pLandau = Vectorize(function(x,mu=log(pi/2),sigma=pi/2,log=FALSE,lower.tail=TRUE) {
    param = setParam(alpha=1, location=mu, logscale=log(sigma), pm=0)
    return(pEstable(x,param,log=log,lower.tail=lower.tail))
})
qLandau = Vectorize(function(p,mu=log(pi/2),sigma=pi/2,log=FALSE,lower.tail=TRUE) {
    param = setParam(alpha=1, location=mu, logscale=log(sigma), pm=0)
    return(qEstable(p,param,log=log,lower.tail=lower.tail))
})
rLandau = Vectorize(function(n,mu=log(pi/2),sigma=pi/2) {
    return(qLandau(runif(n),mu,sigma))
})
dharmonicmeanp = Vectorize(function(x, L, log=FALSE) {
    x=pmax(1e-300,x); # Would be better to calculate limit
    if(log) return(dLandau(1/x, mu=log(L)+1+psigamma(1)-log(2/pi), sigma=pi/2, log=TRUE)-2*log(x))
    return(dLandau(1/x, mu=log(L)+1+psigamma(1)-log(2/pi), sigma=pi/2, log=FALSE)/x^2)
})
pharmonicmeanp = Vectorize(function(x, L, log=FALSE, lower.tail=TRUE) {
    return(pLandau(1/x, mu=log(L)+1+psigamma(1)-log(2/pi), sigma=pi/2, log=log, lower.tail=!lower.tail))
})
qharmonicmeanp = Vectorize(function(p, L, log=FALSE, lower.tail=TRUE) {
    return(1/qLandau(p, mu=log(L)+1+psigamma(1)-log(2/pi), sigma=pi/2, log=log, lower.tail=!lower.tail))
})
rharmonicmeanp = Vectorize(function(n, L) {
    return(qharmonicmeanp(runif(n),L))
})
dmamml = Vectorize(function(x, L, df, log=FALSE) {
    c = ifelse(df==2,1,ifelse(df<2,x*(1-pgamma(log(x),df/2,1)),L*x*(1-pgamma(log(L*x),df/2,1))))
    return(dLandau(x, mu=c*(log(L)+1+psigamma(1)-log(2/pi)), sigma=c*pi/2, log=log))
})
pmamml = Vectorize(function(x, L, df, log=FALSE, lower.tail=TRUE) {
    c = ifelse(df==2,1,ifelse(df<2,x*(1-pgamma(log(x),df/2,1)),L*x*(1-pgamma(log(L*x),df/2,1))))
    return(pLandau(x, mu=c*(log(L)+1+psigamma(1)-log(2/pi)), sigma=c*pi/2, log=log, lower.tail=lower.tail))
})
qmamml = Vectorize(function(p, L, df, log=FALSE, lower.tail=TRUE, xmin=1+1e-12, xmax=1e12) {
    f = function(x) pmamml(x,L,df,log=log,lower.tail=lower.tail)-p
    return(uniroot(f,c(xmin,xmax))$root)
})
rmamml = Vectorize(function(n, L, df) {
    return(qmamml(runif(n),L,df))
})

