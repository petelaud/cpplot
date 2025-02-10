#' Jeffreys and other approximate Bayesian confidence intervals for a single
#' binomial or Poisson rate.
#'
#' Generalised approximate Bayesian confidence intervals based on a Beta (for
#' binomial rates) or Gamma (for Poisson rates) conjugate priors. Encompassing
#' the Jeffreys method (with Beta(0.5, 0.5) or Gamma(0.5) respectively), as well
#' as any user-specified prior distribution. Clopper-Pearson method (as
#' quantiles of a Beta distribution as described in Brown et al. 2001) also
#' included by way of a "continuity correction" parameter.
#'
#' @param x Numeric vector of number of events.
#' @param n Numeric vector of sample sizes (for binomial rates) or exposure
#'   times (for Poisson rates).
#' @param ai,bi Numbers defining the Beta prior distribution (default ai = bi =
#'   0.5 for Jeffreys interval). Gamma prior for Poisson rates requires only ai.
#' @param cc Number or logical specifying (amount of) "continuity correction".
#'   cc = 0 (default) gives Jeffreys interval, cc = 0.5 gives the
#'   Clopper-Pearson interval (or Garwood for Poisson). A value between 0 and
#'   0.5 allows a compromise between proximate and conservative coverage.
#' @param level Number specifying confidence level (between 0 and 1, default
#'   0.95).
#' @param distrib Character string indicating distribution assumed for the input
#'   data: "bin" = binomial (default), "poi" = Poisson.
#' @param adj Logical (default TRUE) indicating whether to apply the boundary
#'   adjustment recommended on p108 of Brown et al. (set to FALSE if informative
#'   priors are used)
#' @param ... Other arguments.
#' @importFrom stats qbeta qgamma qnorm
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @examples
#' # Jeffreys method:
#' jeffreysci(x = 5, n = 56)
#' @references
#'   Laud PJ. Equal-tailed confidence intervals for comparison of
#'   rates. Pharmaceutical Statistics 2017; 16:334-348.
#'
#'   Brown LD, Cai TT, DasGupta A. Interval estimation for a binomial
#'   proportion. Statistical Science 2001; 16(2):101-133
#' @export
jeffreysci <- function(x,
                       n,
                       ai = 0.5,
                       bi = 0.5,
                       cc = 0,
                       level = 0.95,
                       distrib = "bin",
                       adj = TRUE,
                       ...) {
  if (!is.numeric(c(x, n))) {
    print("Non-numeric inputs!")
    stop()
  }
  if (any(c(x, n) < 0)) {
    print("Negative inputs!")
    stop()
  }
  if (distrib == "bin" && (any(x > n + 0.001))) {
    print("x > n not possible for distrib = 'bin'")
    stop()
  }
  if (as.character(cc) == "TRUE") cc <- 0.5

  alpha <- 1 - level
  if (distrib == "bin") {
    CI_lower <- qbeta(alpha / 2, x + (ai - cc), n - x + (bi + cc))
    est <- qbeta(0.5, x + (ai), n - x + (bi)) # Obtain phat as the median
    CI_upper <- qbeta(1 - alpha / 2, x + (ai + cc), n - x + (bi - cc))
    if (adj == TRUE) {
      # adjustment at boundary values
      CI_lower[x == 0] <- 0
      CI_upper[x == n] <- 1
      est[x == 0] <- 0
      est[x == n] <- 1
    }
  } else if (distrib == "poi") {
    # Jeffreys prior for Poisson rate uses gamma distribution,
    # as defined in Li et al. with "continuity correction" from Laud 2017.
    CI_lower <- qgamma(alpha / 2, x + (ai - cc), scale = 1 / n)
    est <- qgamma(0.5, x + (ai), scale = 1 / n)
    if (adj == TRUE) {
      CI_lower[x == 0] <- 0
      est[x == 0] <- 0
    }
    CI_upper <- qgamma(1 - alpha / 2, (x + (ai + cc)), scale = 1 / n)
  }
  CI <- cbind(Lower = CI_lower, est = est, Upper = CI_upper)
  CI
}



#' Wilson score interval, and equivalent Rao score for Poisson data
#'
#' with optional continuity correction
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#' Altman DG, Machin D, Bryant TN et al (2000) Statistics with confidence,
#' 2nd edn. BMJ Books, Bristol
#'
#' Li et al 2014. Comput Stat (2014) 29:869–889
#'
#' Labelled as "Second Normal" in REVSTAT – Statistical Journal
#' Volume 10, Number 2, June 2012, 211–227
#' (which provides the continuity correction formula)
#'
#' Schwertman, N.C. and Martinez, R.A. (1994). Approximate Poisson confidence
#' limits, Communication in Statistics — Theory and Methods, 23(5), 1507-1529.
#'
#' @noRd
wilsonci <- function(x,
                     n,
                     level = 0.95,
                     cc = FALSE,
                     distrib = "bin") {
  if (as.character(cc) == "TRUE") cc <- 0.5
  corr <- cc / n
  z <- qnorm(1 - (1 - level) / 2)
  est <- x / n
  if (distrib == "bin") {
    lower <- (2 * (x - cc) + z^2 -
                z * sqrt(z^2 - 2 * (2 * cc + cc / n) +
                           4 * ((x / n) * (n * (1 - x / n) + 2 * cc)))
    ) / (2 * (n + z^2))
    lower[x == 0] <- 0 # See Newcombe 1998
    upper <- (2 * (x + cc) + z^2 +
                z * sqrt(z^2 + 2 * (2 * cc - cc / n) +
                           4 * ((x / n) * (n * (1 - x / n) - 2 * cc)))
    ) / (2 * (n + z^2))
    upper[x == n] <- 1
  } else if (distrib == "poi") {
    lower <- ((x - cc) + z^2 / 2 - z * sqrt(x - cc + z^2 / 4)) / n
    lower[x == 0] <- 0
    upper <- ((x + cc) + z^2 / 2 + z * sqrt(x + cc + z^2 / 4)) / n
  }
  cbind(Lower = lower, MLE = est, Upper = upper)
}


#' Bisection root-finding
#'
#' vectorized limit-finding routine - turns out not to be any quicker but is
#' neater. The bisection method is just as efficient as the secant method
#' suggested by G&N, and affords greater control over whether the final estimate
#' has score<z the secant method is better for RR and for Poisson rates, where
#' there is no upper bound for d, however it is not guaranteed to converge New
#' version not reliant on point estimate This could be modified to solve upper
#' and lower limits simultaneously
#'
#' @inheritParams scoreci
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#'
#' @noRd
bisect <- function(ftn,
                   contrast,
                   distrib,
                   precis,
                   max.iter = 100,
                   uplow = "low") {
  tiny <- (10^-(precis)) / 2
  nstrat <- length(eval(ftn(1)))
  hi <- rep(1, nstrat)
  lo <- rep(-1, nstrat)
  dp <- 2
  niter <- 1
  while (niter <= max.iter && any(dp > tiny | is.na(hi))) {
    dp <- 0.5 * dp
    mid <- pmax(-1, pmin(1, round((hi + lo) / 2, 10)))
    # rounding avoids machine precision problem with, e.g. 7/10-6/10
    if (contrast == "RD" && distrib == "bin") {
      scor <- ftn(mid)
    } else if (contrast == "RD" && distrib == "poi") {
      scor <- ftn(round(tan(pi * mid / 2), 10))
      # avoid machine precision producing values outside [-1, 1]
    } else if (contrast %in% c("RR", "OR") ||
               (contrast == "p" && distrib == "poi")) {
      scor <- ftn(round(tan(pi * (mid + 1) / 4), 10))
      # avoid machine precision producing values outside [-1, 1]
    } else if (contrast == "p" && distrib == "bin") {
      scor <- ftn((mid + 1) / 2)
    }
    check <- (scor <= 0) | is.na(scor)
    # ??scor=NA only happens when |p1-p2|=1 and |theta|=1 for RD
    # (in which case hi==lo anyway), or if p1=p2=0
    # also for RR when p1=0 and theta=0
    hi[check] <- mid[check]
    lo[!check] <- mid[!check]
    niter <- niter + 1
  }
  if (uplow == "low") {
    best <- lo
  } else {
    best <- hi
  }
  if (contrast == "RD" && distrib == "bin") {
    return(best)
  } else if ((contrast %in% c("RD") && distrib == "poi")) {
    return(tan(best * pi / 2))
  } else if (contrast %in% c("RR", "OR") ||
             (contrast == "p" && distrib == "poi")) {
    return(tan((best + 1) * pi / 4))
  } else if (contrast == "p" && distrib == "bin") {
    return((best + 1) / 2)
  }
}

#' Internal function - superseded by new wilsonci function
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#'
#' @noRd
quadroot <- function(a, b, c_) {
  # GET ROOTS OF A QUADRATIC EQUATION
  r1x <- (-b + sqrt(b^2 - 4 * a * c_)) / (2 * a)
  r2x <- (-b - sqrt(b^2 - 4 * a * c_)) / (2 * a)
  r1 <- pmin(r1x, r2x)
  r2 <- pmax(r1x, r2x)
  cbind(r1, r2)
}
