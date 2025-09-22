

# For a given parameter combination: p1, p2, (psi or phi)
# calculate the probability density function for each combination x = c(a, b, c, d)
# vectorized so x can be input as an array with 4 columns
pdfpair <- function(x,
                    p1 = NULL,
                    p2 = NULL,
                    phi = NULL,
                    psi = NULL
) {
  if (!is.null(phi)) {
    p11 <- p1 * p2 + phi * sqrt(p1 * (1 - p1) * p2 * (1 - p2)) # with parameterisation using phi
  }
  # Or using psi=(p11*p22)/(p12*p21) (labelled as theta in Fagerland)
  if (!is.null(psi)) {
    A <- psi - 1
    B <- p1 + p2 - 1 - psi * (p1 + p2)
    C <- psi * p1 * p2
    p11 <- ifelse(psi == 1, -C / B,(-B - sqrt(B^2 - 4 * A * C)) / (2 * A))
  }
  p12 <- p1 - p11
  p21 <- p2 - p11
  p22 <- pmax(0, 1 - p11 - p12 - p21)
  dens <- exp(
    lfactorial(rowSums(x)) -
      (lfactorial(x[, 1]) + lfactorial(x[, 2]) + lfactorial(x[, 3]) + lfactorial(x[, 4]))
  ) * p11^x[, 1] * p12^x[, 2] * p21^x[, 3] * (1 - p11 - p12 - p21)^(x[, 4])
  dens[(p21 <= 0 & p12 > 0) | (p12 <= 0 & p21 > 0)] <- NA
  return(dens)
}




#' Internal function to evaluate the score at a given value of theta
#' Copied from ratesci package
#'
#' For iterative calculations:
#' function to evaluate the score at a given value of theta, given the observed
#' data for paired binomial RD and RR. uses the MLE solution (and notation)
#' given in Fagerland from Tango (1998/1999) & Tang (2003)
#' This function is not vectorised
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#'   Tango T. Equivalence test and confidence interval for the difference
#'   in proportions for the paired-sample design.
#'   Statistics in Medicine 1998; 17:891-908
#'
#'   Nam J-M, Blackwelder WC. Analysis of the ratio of marginal
#'   probabilities in a matched-pair setting.
#'   Stat Med 2002; 21(5):689–699
#'
#'   Tang N-S, Tang M-L, Chan ISF. On tests of equivalence via non-unity
#'   relative risk for matched-pair design.
#'   Statistics in Medicine 2003; 22:1217-1233
#'
#' @inheritParams pairbinci
#'
#' @noRd
scorepair <- function(theta,
                      x,
                      contrast = "RD",
                      cc = FALSE,
                      cctype = "new",
                      skew = FALSE,
                      bcf = FALSE,
                      ...) {
  N <- sum(x)
  lambda <- switch(as.character(bcf),
                   "TRUE" = N / (N - 1),
                   "FALSE" = 1
  )
  if (as.character(cc) == "TRUE") cc <- 0.5

  if (contrast == "RD") {
    # notation per Tango 1999 letter, divided by N
    # (to get numerator on the right scale for skewness correction)
    # Variance is divided by N^2 accordingly
    # and continuity correction term also divided by N
    Stheta <- ((x[2] - x[3]) - N * theta) / N
    A <- 2 * N
    B <- -x[2] - x[3] + (2 * N - x[2] + x[3]) * theta
    C_ <- -x[3] * theta * (1 - theta)
    num <- (-B + Re(sqrt(as.complex(B^2 - 4 * A * C_))))
    p21 <- ifelse(num == 0, 0, num / (2 * A))
    corr <- 2 * cc * sign(Stheta) / N
    p12 <- p21 + theta
    # From Tango
    p11 <- ifelse(x[1] == 0, 0, x[1] / (x[1] + x[4]) * (1 - p12 - p21))
    p22 <- 1 - p11 - p12 - p21
    p2d <- pmin(1, pmax(0, p21 + p11))
    p1d <- p2d + theta

    # Tango variance
    #    V <- pmax(0, N * (2 * p21 + theta * (1 - theta))) * lambda / (N^2)
    # Equivalent
    V <- pmax(0, (p1d * (1 - p1d) + p2d * (1 - p2d) -
                    2 * (p11 * p22 - p12 * p21)) / N) * lambda
    mu3 <- (p1d * (1 - p1d) * (1 - 2 * p1d) +
              ((-1)^3) * p2d * (1 - p2d) * (1 - 2 * p2d) +
              3 * (-1) * (p11 * (1 - p1d)^2 + p21 * p1d^2 - p1d * p2d * (1 - p1d)) +
              3 * ((-1)^2) * (p11 * (1 - p2d)^2 + p12 * p2d^2 - p1d * p2d * (1 - p2d))
    ) / (N^2)
  }
  if (contrast == "RR") {
    # per Tang 2003, but divided by N
    Stheta <- ((x[2] + x[1]) - (x[3] + x[1]) * theta) / N
    A <- N * (1 + theta)
    B <- (x[1] + x[3]) * theta^2 - (x[1] + x[2] + 2 * x[3])
    C_ <- x[3] * (1 - theta) * (x[1] + x[2] + x[3]) / N
    num <- (-B + Re(sqrt(as.complex(B^2 - 4 * A * C_))))
    q21 <- ifelse(num == 0, 0, num / (2 * A))

    if (cctype == "constant") corr <- cc * 2 * sign(Stheta) / N
    if (cctype == "delrocco") corr <- cc * (x[1] + x[3]) / N * sign(Stheta) / N
    # Equivariant continuity correction for RR, aligned with McNemar cc.
    if (cctype == "new") corr <- cc * (1 + theta) * sign(Stheta) / N

    q12 <- (q21 + (theta - 1) * (1 - x[4] / N)) / theta
    # Below from Tang 2003
    q11 <- (1 - x[4] / N - (1 + theta) * q21) / theta
    q22 <- 1 - q11 - q12 - q21
    p2d <- q21 + q11
    #    p1d <- q12 + q11
    p1d <- p2d * theta

    # Tang variance
    #    V <- pmax(0, N * (1 + theta) * q21 + (x[1] + x[2] + x[3]) * (theta - 1)) * lambda / (N^2)
    # Nam-Blackwelder variance
    #    V <- pmax(0, theta*(q12 + q21) / N) * lambda
    # Equivalent consistent with scoreci notation
    V <- pmax(0, (p1d * (1 - p1d) + theta^2 * p2d * (1 - p2d) -
                    2 * theta * (q11 * q22 - q12 * q21)) / N) * lambda
    mu3 <- (p1d * (1 - p1d) * (1 - 2 * p1d) +
              ((-theta)^3) * p2d * (1 - p2d) * (1 - 2 * p2d) +
              3 * (-theta) * (q11 * (1 - p1d)^2 + q21 * p1d^2 - p1d * p2d * (1 - p1d)) +
              3 * ((-theta)^2) * (q11 * (1 - p2d)^2 + q12 * p2d^2 - p1d * p2d * (1 - p2d))) / (N^2)
  }
  scterm <- mu3 / (6 * V^(3 / 2))
  scterm[abs(mu3) < 1E-10] <- 0 # Avoids issues with e.g. x = c(1, 3, 0, 6)
  score1 <- ifelse(Stheta == 0, 0, (Stheta - corr) / sqrt(V))
  A <- scterm
  B <- 1
  C_ <- -(score1 + scterm)
  num <- (-B + sqrt(pmax(0, B^2 - 4 * A * C_)))
  score <- ifelse((skew == FALSE | scterm == 0),
                  score1, num / (2 * A)
  )
  score[abs(Stheta) < abs(corr)] <- 0

  pval <- pnorm(score)
  outlist <- list(
    score = score, p1d = p1d, p2d = p2d,
    mu3 = mu3, pval = pval
  )
  return(outlist)
}



