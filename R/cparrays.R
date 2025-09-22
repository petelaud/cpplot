# Evaluation of confidence intervals for RD or RR from paired data
# 2025 update using phi correlation parameter
# and adding MOVER methods

rm(list = ls())
#install.packages("zoo")
pak::pak('petelaud/ratesci')

set.seed(2012) #ensure we use the same jitters for each run

root <- "/Users/ssu/Documents/"
outpath <- paste(root, "Main/Courses_papers/skewscore/paired/", sep = "")

mysum <- function(x) sum(x,na.rm=T)
mymean <- function(x) mean(x,na.rm=T)

# Conversion of parameters for RR
# convert phi to psi, or vice versa, for a given p1 and (p2 or theta)
# also output p11, p12, p21, p22 for a given p1, p2 and (phi or psi)
params <- function(p1,
                   p2 = NULL,
                   theta = NULL,
                   contrast = "RD",
                   psi = NULL,
                   phi = NULL
                   ) {
  if (is.null(p2)) {
    if (contrast == "RR") p2 <- p1 / theta
    if (contrast == "RD") p2 <- p1 - theta
  }
  if (is.null(phi)) {
    A <- psi-1
    B <- p1+p2-1-psi*(p1+p2)
    C <- psi*p1*p2
    p11 <- (-B - sqrt(B^2-4*A*C))/(2*A)
  }
  if (is.null(psi)) {
    p11 <- p1 * p2 + phi*sqrt(p1*(1-p1)*p2*(1-p2))
  }
  p12 <- p1 - p11
  p21 <- p2 - p11
  p22 <- 1 - p11 - p12 - p21

  if (is.null(phi)) {
    phi <- (p11-p1*p2)/((p1*(1-p1)*p2*(1-p2)))^0.5
  }
  if (is.null(psi)) {
    psi <- p11*p22/(p12*p21)
  }
  cbind(p1, p2, psi, phi, p11, p12, p21, p22)
}

# Plot to illustrate relationship between phi and psi
if (FALSE) {
p1s <- seq(0,1,0.05)

dev.off()
# Show that Fagerland's selected psi values translate to low correlation of 0.1-0.3
# Plot phi vs p1, with RD=0.2 and for psi = 2 or 4
plot(p1s, params(p1=p1s, theta=0, psi=2, contrast="RD")[,4], ylim = c(0, 0.9), type='n',
     main=expression(paste(theta["RD"], " = 0.2")), xlab = "p1",
     ylab = expression(phi))
lines(p1s, params(p1=p1s, theta=0.2, psi=2)[,4], lty=2)
lines(p1s, params(p1=p1s, theta=0.2, psi=4)[,4], lty=1)
abline(h=0.2, lty=3)
# Plot phi vs p1, with RR=4 and for psi = 2 or 4
plot(p1s, params(p1=p1s, theta=1, psi=2, contrast="RR")[,4], ylim = c(0, 0.9), type='n',
     main=expression(paste(theta["RR"], " = 4")), xlab = "p1",
     ylab = expression(phi))
lines(p1s, params(p1=p1s, theta=4, psi=2, contrast="RR")[,4], lty=2)
lines(p1s, params(p1=p1s, theta=4, psi=4, contrast="RR")[,4], lty=1)
abline(h=0.2, lty=3)
#lines(p1s, params(p1=p1s, theta=0, psi=10)[,4], lty=2)
#lines(p1s, params(p1=p1s, theta=0, psi=100)[,4], lty=3)
del <- 0.1
plot(p1s[p1s >= del], params(p1s[p1s >= del], theta=0)[,4], ylim = c(0, 0.9), type='n',
     main=expression(paste(theta["RD"], " = 0.1")), xlab = "p1",
     ylab = expression(phi))
lines(p1s, params(p1=p1s, theta=del, psi=2)[,4], lty=2)
lines(p1s, params(p1=p1s, theta=del, psi=4)[,4], lty=1)
del <- 0.2
plot(p1s[p1s >= del], params(p1s[p1s >= del], theta=0)[,4], ylim = c(0, 0.9), type='n',
     main=expression(paste(theta["RD"], " = 0.2")), xlab = "p1",
     ylab = expression(phi))
lines(p1s, params(p1=p1s, theta=del, psi=2)[,4], lty=2)
lines(p1s, params(p1=p1s, theta=del, psi=4)[,4], lty=1)
max(params(p1=p1s, theta=del, psi=2)[,4], na.rm=T)
max(params(p1=p1s, theta=0, psi=4)[,4], na.rm=T)

# Original version - triple plot for psi=2,10,100 and RR=1,2,4
par(mfrow=c(1,3), pty='s', mar=(c(2,4,2,0)+0.1))
plot(p1s, params(p1=p1s, theta=1, psi=2)[,4], ylim = c(0, 0.9), type='n',
     main=expression(paste(theta["RR"], " = 1")), xlab = "p1",
     ylab = expression(phi))
lines(p1s, params(p1=p1s, theta=1, psi=2)[,4], lty=1)
lines(p1s, params(p1=p1s, theta=1, psi=4)[,4], lty=1)
lines(p1s, params(p1=p1s, theta=1, psi=10)[,4], lty=2)
lines(p1s, params(p1=p1s, theta=1, psi=100)[,4], lty=3)
plot(p1s, params(p1=p1s, theta=4, psi=2)[,4], ylim = c(0, 0.9), type='n',
     main=expression(paste(theta["RR"], " = 2")), xlab = "p1",
     ylab = expression(phi))
lines(p1s, params(p1=p1s, theta=2, psi=2)[,4], lty=1)
lines(p1s, params(p1=p1s, theta=2, psi=10)[,4], lty=2)
lines(p1s, params(p1=p1s, theta=2, psi=100)[,4], lty=3)
plot(p1s, params(p1=p1s, theta=4, psi=2)[,4], ylim = c(0, 0.9), type='n',
     main=expression(paste(theta["RR"], " = 4")), xlab = "p1",
     ylab = expression(phi))
lines(p1s, params(p1=p1s, theta=4, psi=2)[,4], lty=1)
lines(p1s, params(p1=p1s, theta=4, psi=10)[,4], lty=2)
lines(p1s, params(p1=p1s, theta=4, psi=100)[,4], lty=3)
mtext(expression(paste("Relationship between ", phi, " and ", psi, " for ", psi,
                       " = 2, 10, 100")),
      side=3, line=-4, outer=TRUE, cex=1.5)
mtext(expression(paste("at selected values of ", theta["RR"])),
      side=3, line=-6, outer=TRUE, cex=1.5)
}


# Calculate all (or selected) CI methods for a given set of x's
allpairci <- function(xs,
                      contrast = "RD",
                      alpha = 0.05,
                      methods = NULL
                      ) {
  # Force vector input to an array
  if (is.vector(xs)) dim(xs) <- c(1, length(xs))
  lenxs <- dim(xs)[1]
  if (!is.null(methods)) {
    mymethods <- methods
  } else {
    if (contrast %in% c("RD", "RR")) {
      mymethods <- c("AS", "AS-bc", "SCAS", "SCAS-bc",
                     "MOVER-W", "MOVER-NW", "MOVER-NJ", "MOVER-NS",
                     "SCAS-c5", "SCAS-c25", "SCAS-c125",
                     "MOVER-c5", "MOVER-c25", "MOVER-c125", "BP")
      if (contrast == "RR") {
        mymethods <- c(mymethods, "BP-J") #, "Tang-ccdr")
      }
    } else if (contrast == "OR") {
      mymethods <- c("SCASp", "SCASpu", "Jeffreys", "mid-p", "Wilson",
                     "SCASp-c5", "SCASp-c25", "SCASp-c125", "C-P",
                     "Jeffreys-c25", "Jeffreys-c125", "Wilson-c", "midp-c25")
    }
  }
  ci <- array(NA, dim = c(lenxs, 2, length(mymethods)))
  dimnames(ci)[[3]] <- mymethods
  if (contrast %in% c("RD", "RR")) {
    # "AS" is Tango asymptotic score method for RD, and Tang for RR
    if ("AS" %in% mymethods) ci[, 1:2, "AS"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "Score_closed", skew = FALSE, bcf = FALSE, level = 1-alpha)$estimates[,c(1,3)]))
    if ("AS-bc" %in% mymethods)  ci[, 1:2, "AS-bc"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "Score_closed", skew = FALSE, bcf=TRUE, level = 1-alpha)$estimates[,c(1,3)]))
    # "SCAS" adds the skewness correction
    if ("SCAS" %in% mymethods)  ci[, 1:2, "SCAS"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "Score", skew=TRUE, bcf = FALSE, level = 1-alpha)$estimates[,c(1,3)]))
    # plus the 'N-1' variance bias correction
    if ("SCAS-bc" %in% mymethods) ci[, 1:2, "SCAS-bc"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "Score", skew=TRUE, bcf=TRUE, level = 1-alpha)$estimates[,c(1,3)]))
    # Explore continuity adjustments
    if ("SCAS-c5" %in% mymethods)  ci[, 1:2, "SCAS-c5"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "Score", skew=TRUE, bcf=TRUE, cc=0.5, level = 1-alpha)$estimates[,c(1,3)]))
    if ("SCAS-c25" %in% mymethods)  ci[, 1:2, "SCAS-c25"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "Score",  skew=TRUE, bcf=TRUE, cc=0.25, level = 1-alpha)$estimates[,c(1,3)]))
    if ("SCAS-c125" %in% mymethods)  ci[, 1:2, "SCAS-c125"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "Score", skew=TRUE, bcf=TRUE, cc=0.125, level = 1-alpha)$estimates[,c(1,3)]))
    # MOVER methods, first without modification
    if ("MOVER-W" %in% mymethods) ci[, 1:2, "MOVER-W"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "MOVER", moverbase = "wilson", level = 1-alpha)$estimates[,c(1,3)]))
    # "-N" adds Newcombe's correlation correction
    if ("MOVER-NW" %in% mymethods)  ci[, 1:2, "MOVER-NW"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "MOVER_newc", moverbase = "wilson", level = 1-alpha)$estimates[,c(1,3)]))
    # then try different input methods: Jeffreys & SCASp
    if ("MOVER-NJ" %in% mymethods) ci[, 1:2, "MOVER-NJ"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "MOVER_newc", moverbase = "jeff", level = 1-alpha)$estimates[,c(1,3)]))
    if ("MOVER-NS" %in% mymethods)  ci[, 1:2, "MOVER-NS"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "MOVER_newc", moverbase = "SCASp", level = 1-alpha)$estimates[,c(1,3)]))
    # Explore continuity adjustments
    if ("MOVER-c5" %in% mymethods)  ci[, 1:2, "MOVER-c5"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "MOVER_newc", moverbase = "jeff", level = 1-alpha, cc=0.5)$estimates[,c(1,3)]))
    if ("MOVER-c25" %in% mymethods)  ci[, 1:2, "MOVER-c25"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "MOVER_newc", moverbase = "jeff", level = 1-alpha, cc=0.25)$estimates[,c(1,3)]))
    if ("MOVER-c125" %in% mymethods)  ci[, 1:2, "MOVER-c125"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "MOVER_newc", moverbase = "jeff", level = 1-alpha, cc=0.125)$estimates[,c(1,3)]))
    # Bonett-Price hybrid method for RR, or their adjusted Wald method for RD
    if ("BP" %in% mymethods) ci[, 1:2, "BP"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "BP", moverbase = "wilson", level = 1-alpha)$estimates[,c(1,3)]))
    if (alpha == 0.05) {
      # Further evaluation of TDAS method proposed for paired analysis in Laud2017
      # doesn't compare favourably with new paired SCAS method
#      if ("TDAS" %in% mymethods)  ci[, 1:2, "TDAS"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "TDAS")$estimates[c(1,3)]))
      # Also stratified SCAS - not good, ABANDON
      # if ("SCASstrat" %in% mymethods)  ci[, 1:2, "SCASstrat"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "SCASstrat")$estimates[c(1,3)]))
    }
    if (contrast == "RR") {
      # DelRocco's version of continuity correction for RR - not equivariant so ruled out
#      if ("Tang-ccdr" %in% mymethods)   ci[, 1:2, "Tang-ccdr"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, level = 1-alpha, cc=0.5, cctype="delrocco")$estimates[,c(1,3)]))
      # Proposed variation of Bonett-Price hybrid method, incorporating Jeffreys intervals
      if ("BP-J" %in% mymethods)    ci[, 1:2, "BP-J"] <- t(sapply(1:lenxs,function(i) {
        if (sum(xs[i, 2:3]) == 0 ) outdata <- c(0, Inf)
        else outdata <- pairbinci(x = xs[i,], contrast = contrast, method = "BP", moverbase = "jeff", level = 1-alpha)$estimates[,c(1,3)]
        outdata
      }
            ))
    }
  } else if (contrast == "OR") {
    # Explore various options for transformed binomial intervals for conditional OR
    # Transformed SCASp with experimental bcf using N/(N-1) to match 'N-1' test for association
    if ("SCASp" %in% mymethods) ci[, 1:2, "SCASp"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "SCASp", bcf = TRUE, level = 1-alpha)$estimates[,c(1,3)]))
    # Transformed Uncorrected SCAS (i.e. skewness-corrected Wilson)
    if ("SCASpu" %in% mymethods) ci[, 1:2, "SCASpu"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "SCASp", bcf = FALSE, level = 1-alpha)$estimates[,c(1,3)]))
    # Transformed Jeffreys equal-tailed interval
    if ("Jeffreys" %in% mymethods) ci[, 1:2, "Jeffreys"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "jeff", level = 1-alpha)$estimates[,c(1,3)]))
    # Transformed Clopper-Pearson mid-p
    if ("mid-p" %in% mymethods) ci[, 1:2, "mid-p"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "midp", level = 1-alpha)$estimates[,c(1,3)]))
    # Transformed Wilson score
    if ("Wilson" %in% mymethods) ci[, 1:2, "Wilson"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "wilson", level = 1-alpha)$estimates[,c(1,3)]))
    # Explore continuity adjustments - could also try reduced gamma variations, e.g. cc=0.125
    if ("C-P" %in% mymethods) ci[, 1:2, "C-P"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "jeff", level = 1-alpha, cc = TRUE)$estimates[,c(1,3)]))
    if ("Wilson-c" %in% mymethods) ci[, 1:2, "Wilson-c"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "wilson", level = 1-alpha, cc = TRUE)$estimates[,c(1,3)]))
    if ("SCASp-c5" %in% mymethods) ci[, 1:2, "SCASp-c5"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "SCASp", bcf = TRUE, cc = 0.5, level = 1-alpha)$estimates[,c(1,3)]))
    if ("SCASp-c25" %in% mymethods) ci[, 1:2, "SCASp-c25"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "SCASp", bcf = TRUE, cc = 0.25, level = 1-alpha)$estimates[,c(1,3)]))
    if ("SCASp-c125" %in% mymethods) ci[, 1:2, "SCASp-c125"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "SCASp", bcf = TRUE, cc = 0.125, level = 1-alpha)$estimates[,c(1,3)]))
    # Transformed continuity-adjusted Jeffreys equal-tailed interval
    if ("Jeffreys-c5" %in% mymethods) ci[, 1:2, "Jeffreys-c5"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "jeff", cc = 0.5, level = 1-alpha)$estimates[,c(1,3)]))
    if ("Jeffreys-c25" %in% mymethods) ci[, 1:2, "Jeffreys-c25"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "jeff", cc = 0.25, level = 1-alpha)$estimates[,c(1,3)]))
    if ("Jeffreys-c125" %in% mymethods) ci[, 1:2, "Jeffreys-c125"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "jeff", cc = 0.125, level = 1-alpha)$estimates[,c(1,3)]))
    # Transformed continuity-adjusted midp interval
    if ("midp-c25" %in% mymethods) ci[, 1:2, "midp-c25"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method = "midp", cc = 0.25, level = 1-alpha)$estimates[,c(1,3)]))
    # (Add conditional logistic regression CI?)
  }
  dimnames(ci)[[2]] <- c("LCL", "UCL")
  ci
}

# Calculates all possible confidence intervals for a given n,
# save the results to a data file
cifun <- function( n,
                   contrast = "RD",
                   alph = c(0.01, 0.05, 0.1),
                   methods = NULL,
                   prerun = FALSE
                   ) {

  g <- (expand.grid(x11 = 0:n, x12 = 0:n, x21 = 0:n))
  # reduce to possible combinations of a,b,c for paired data.
  g <- g[(g$x12 <= n - g$x11) &
           (g$x21 <= n - g$x11 - g$x12), ]
  xs <- data.matrix(cbind(g, x22 = n - rowSums(g)))
  lenxs <- dim(xs)[1]
  row.names(xs) <- NULL
  # (Could achieve the same effect by expanding a,b,c,d and selecting all rows where total=n.)
  # 	g <- expand.grid(a=0:n,b=0:n,c=0:n,d=0:n)
  # 	g <- g[rowSums(g) == n,1:3]

  if (!is.null(methods)) {
    mymethods <- methods
  } else {
    tester <- allpairci(x = rep(10,4), contrast = contrast)
    mymethods <- dimnames(tester)[[3]]
  }
  nmeth <- length(mymethods)

  cis <- array(NA, dim = c(lenxs, 2, nmeth, length(alph), 1, 1))
  dimnames(cis)[[2]] <- c("LCL", "UCL")
  dimnames(cis)[[3]] <- mymethods
  dimnames(cis)[[4]] <- paste(100*(1-alph))
  dimnames(cis)[[5]] <- paste(n)
  dimnames(cis)[[6]] <- contrast

  cat("Step 1 of 3: Calculate selected CIs for all possible outcomes for the selected N\n")
  for (alpha in alph) {
    cat(paste0("alpha=", alpha,"\n"))
    # For the given sample size N, we calculate confidence intervals,
    # for each method, for every possible combination of observed frequencies
#    ci <- allpairci(xs = xs, contrast = contrast, alpha = alpha, methods = mymethods)
#    cis[, , , paste(100 * (1 - alpha)), 1, 1] <- ci
    # Process one x at a time so we can show progress and expected time
    pb <- pbapply::timerProgressBar(min = 0, max = dim(xs)[1], width = 50, char = '=', style = 1) #timer bar
    on.exit(pbapply::closepb(pb))
    for (i in 1:dim(xs)[[1]]) {
      pbapply::setTimerProgressBar(pb, i)
      ci <- allpairci(xs = xs[i, ], contrast = contrast, alpha = alpha, methods = mymethods)
      cis[i, , , paste(100 * (1 - alpha)), 1, 1] <- ci
    }
    pbapply::closepb(pb)
  }

  ciarrays <- list(xs = xs, cis = cis)
  save(ciarrays, file = paste0(outpath, "cis.", contrast, ".", n, ".Rdata"))
  ciarrays
}

# Take array of CIs from cifun(), create grid of PSPs and calculate
# coverage probabilities (raw and smoothed) and summaries
cpfun <- function(
                  ciarrays = myciarrays,
                  alph = NULL,
                  psis = NULL,
                  phis = c(0, 0.25, 0.5, 0.75),
                  n.grid = 200,
                  res.factor = 6,
                  methods = "All",
                  limits = c(0, 1),
                  window = 0.08,
                  square = TRUE,
                  smooth = TRUE,
                  prerun = F,
                  jitt = TRUE,
                  sided = "R",
                  pcut = 1E-10
                  ) {

  xs <- ciarrays[["xs"]]
  cis <- ciarrays[["cis"]]
  contrast <- dimnames(cis)[[6]]
  mymethods <- longlab <- dimnames(cis)[[3]]
  if (is.null(alph)) alph <- 1 - (as.numeric(dimnames(cis)[[4]])/100)
  nmeth <- length(mymethods)
  n <- as.numeric(dimnames(cis)[[5]])

  if (!is.null(psis)) par3 <- psis
  if (!is.null(phis)) par3 <- phis
  p.grid <- limits[1] + ppoints(n.grid, a = 0.5) * (limits[2] - limits[1])
  nn <- qbinom(0.99999, n, limits[2]) # Upper limit for n when considering Poisson rates?
  p1 <- p2 <- p.grid

  mastercp <- array(NA, dim = c(length(p1), length(p2), length(par3),
                                length(mymethods), length(alph), 10, 1, 1))
  dimnames(mastercp) <- list(paste(p1),
                             paste(p2),
                             paste(par3),
                             mymethods,
                             paste(100*(1-alph)),
                             c("cp", "lncp", "rncp", "dncp", "avecp",
                               "avelncp", "averncp", "len", "locindex", "avelocindex"),
                             paste(n),
                             contrast)

  summarylist <- c("meanCP", "minCP", "pctCons", "pctBad", "pctBad2", "pctnear",
                   "pctAveCons", "pctAvenear", "pctAvenearlo", "pctAvenearhi", "pctAveBad",
                   "meanMNCP", "meanDNCP", "maxLNCP", "pctCons.1side", "pctBad.1side", "pctBad.DNCP",
                   "pctnear.1side", "pctAvenear.1side", "typeI", "maxtypeI", "avetypeI",
                   "meanlen", "meanlocindex", "pctgoodloc",
                   "pctnear.DNCP", "pctCons.DNCP", "pctBad.either", "pctCons.both")
  summaries <- array(NA, dim = c(length(par3), length(mymethods), length(alph), length(summarylist), 1, 1))
  dimnames(summaries) <-
    list(paste(par3),
         mymethods,
         paste(100*(1-alph)),
         summarylist,
         paste(n),
         contrast)

  arrays <- list(xs = xs, cis = cis, mastercp = mastercp, summaries = summaries)

  n.grid <- length(p1)
#  set.seed(2012) # ensure we use the same jitters for each run
  px <- pxg <- expand.grid(p1, p2, par3)

  refine <- 1 / n.grid
  if (jitt) {
    jitter <- (limits[2] - limits[1]) * array(runif(2 * dim(pxg)[1], -refine / 2, refine / 2), dim = dim(pxg)) # OK if n.grid is even
    px[, 1:2] <- pxg[, 1:2] + jitter[, 1:2]
  } else {
    px <- pxg
  }

  # theta estimate for every p1,p2 pair
  if (contrast %in% c("RR")) {
    theta <- (px[, 1] / px[, 2])
  } else if (contrast == "OR") {
    # Note this is the conditional OR, p12/p21, not the marginal OR
    # So need to get p12 and p21 from the PSP values of p1, p2 and phi/psi
    if (!is.null(phis)) pij <- params(p1 = px[, 1], p2 = px[, 2], phi = px[, 3])[, 5:8]
    if (!is.null(psis)) pij <- params(p1 = px[, 1], p2 = px[, 2], psi = px[, 3])[, 5:8]
    theta <- pij[, 2] / pij[, 3] # Conditional OR
  } else if (contrast == "RD") {
    theta <- px[, 1] - px[, 2]
  }

  # can we improve efficiency here? -
  # - not without creating an array which is too large for available RAM.

  cat(paste0("Step 2 of 3: calculate coverage probabilities for each alpha\n"))

  for (alpha in alph) {
    cat(paste0("alpha=", alpha,"\n"))
    cat(paste0("Coverage probabilities\n"))
    ci <- cis[,,, paste(100 * (1 - alpha)),,]

    # To obtain coverage probability for n,p1,p2, we determine the probability
    # of each possible pair of observed frequencies a,b given p1,p2 using
    # binomial probabilities
    # and sum over all combinations where p1-p2 is included in the corresponding
    #  CI, for each CI method.
    cpl <- lncpl <- rncpl <- lenl <- locindexl <-
      array(NA, dim = c(dim(px)[1], nmeth))

    # NB method using sapply was no quicker than this loop:
    # NB this is the most time-consuming part of the process when n is small
    # i <- 40 #gives negative probs at p1=0.95, p2=0.35. phi=0.5
    # - pdfpair function updated to give prob=NA for such impossible combinations
    # First switched to parameterisation using psi instead.
    # Then revisited in order to produce plots based on phi (more interpretable)
    pb <- pbapply::timerProgressBar(min = 0, max = dim(px)[1], width = 50, char = '=', style = 1) #timer bar
    on.exit(pbapply::closepb(pb))
    for (i in 1:dim(px)[1]) {
      pbapply::setTimerProgressBar(pb, i)

      if (!is.null(psis)) {
        prob <- pdfpair(p1 = px[i, 1],
                         p2 = px[i, 2],
                         psi = px[i, 3],
                         x = xs)
      } else if (!is.null(phis)) {
        prob <- pdfpair(p1 = px[i, 1],
                        p2 = px[i, 2],
                        phi = px[i, 3],
                        x = xs)
      }

      # can we improve efficiency here? -
      # not without creating an array which is too large for available RAM.

      ########################################
      # - could loop by px and then by alpha?
      # would need to increase dimensions of cpl etc
      ########################################

tryCatch(
  error = function(cnd) {
      message(paste("An error occurred for item", i, px[i, ],":\n"), cnd)

  },
      if (any(!is.na(prob))) {
        cisub <- ci[prob > pcut, , ]
        probsub <- prob[prob > pcut]
        cpl[i, ] <- t(cisub[, 1, ] <= theta[i] & cisub[, 2, ] >= theta[i] & cisub[, 2, ] > cisub[, 1, ]) %*% probsub # 2-sided coverage probability. NB degenerate intervals excluded
        lncpl[i, ] <- t(cisub[, 1, ] > theta[i] | cisub[, 2, ] == cisub[, 1, ]) %*% probsub # L-sided non-coverage (R-side is a mirror image)
        rncpl[i, ] <- t(cisub[, 2, ] < theta[i] | cisub[, 2, ] == cisub[, 1, ]) %*% probsub # R-sided non-coverage

        # not really interested in length, I think location is more important
        # - and its complicated to explain on log scale
        # but included as reviewers might request it
        lens <- cisub[, 2, ] - cisub[, 1, ]
        ## you get whacky results with length on linear scale
        # 		if(contrast %in% c("RR")) lens<-(ci[,2,])-(ci[,1,])
        # 		lens[lens>100]<-100 #workaround for infinite lengths with RR
        if (contrast %in% c("RR")) {
          # for RR, use length on the log scale. Still an issue with infinite lengths though
          lens <- log(cisub[, 2, ]) - log(pmax(0.0000000001, cisub[, 1, ]))
          # Try Newcombe's book suggestion to use U/(1+U) - L/(1+L) - also doesn't resolve the issue
          #lens <- ci[, 2, ]/(1 + ci[, 2, ]) - ci[, 1, ]/(1 + ci[, 1, ])
        }
        lens[lens > 10] <- 10 # workaround for infinite lengths with RR
        lens[cisub[, 2, ] == cisub[, 1, ]] <- 0
        lenl[i, ] <- t(lens) %*% probsub
      }
  ) # End trycatch bracket
    }
    pbapply::closepb(pb)

    mncpl <- rncpl
    thetagt <- px[,1] > px[,2]
    mncpl[thetagt, ] <- lncpl[thetagt, ]

    # Free up some memory
    rm(ci)
    rm(cisub)

    # convert CPs etc to a square grid of p1,p2 points
    mydim <- c(length(p1), length(p2), length(par3), nmeth)
    cp <- array(cpl, dim = mydim)
    lncp <- array(lncpl, dim = mydim)
    rncp <- array(rncpl, dim = mydim)
    mncp <- array(mncpl, dim = mydim)
    dncp <- 1 - cp - mncp
    len <- array(lenl, dim = mydim)
    # Calculate Newcombe's location index
    locindex <- ifelse(cp == 1, 0, mncp / (1 - cp))
    dimnames(cp) <- dimnames(locindex) <- dimnames(lncp) <-
      dimnames(rncp) <- dimnames(len) <- dimnames(dncp) <- dimnames(mncp) <-
      list(paste(p1), paste(p2), paste(par3), mymethods)

    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "cp", , ] <- cp
    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "lncp", , ] <- lncp
    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "rncp", , ] <- rncp
    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "dncp", , ] <- dncp
    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "locindex", , ] <- locindex
    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "len", , ] <- len
    maxncp <- pmax(lncp, rncp)

    meanlen <- format(round(apply(len, 3:4, mymean), 3), nsmall = 3)
    meanCP <- format(round(apply(cp, 3:4, mymean), 3), nsmall = 3)
    minCP <- format(round(apply(cp, 3:4, function(x) min(x, na.rm = TRUE)), 3), nsmall = 3)
    pctCons <- format(100 * round(apply(cp, 3:4, FUN = function(x) mymean(x > (1 - alpha))), 3), nsmall = 1)
    pctBad <- format(100 * round(apply(cp, 3:4, FUN = function(x) mymean(x < (1 - (alpha * 1.1)))), 3), nsmall = 1)
    pctBad2 <- format(100 * round(apply(cp, 3:4, FUN = function(x) mymean(x < (1 - (alpha * 1.2)))), 3), nsmall = 1)
    pctnear <- format(100 * round(apply(cp, 3:4, FUN = function(x) mymean(x > (1 - 1.1 * alpha) & x < (1 - 0.9 * alpha))), 3), nsmall = 1)
    meanlocindex <- format(round(apply(locindex, 3:4, mymean), 3), nsmall = 3)
    pctgoodloc <- format(100 * round(apply(locindex, 3:4, FUN = function(x) mymean(x > 0.4 & x < 0.6)), 3), nsmall = 1)
    # one-sided summaries
    meanMNCP <- format(round(apply(mncp, 3:4, mymean), 3), nsmall = 3)
    meanDNCP <- format(round(apply((1 - cp - mncp), 3:4, mymean), 3), nsmall = 3)
    maxLNCP <- format(round(apply(lncp, 3:4, function(x) max(x, na.rm = TRUE)), 3), nsmall = 3)
    pctBad.1side <- format(100 * round(apply(lncp, 3:4, FUN = function(x) mymean(x > (1.2 * alpha / 2))), 3), nsmall = 1)
    pctBad.DNCP <- format(100 * round(apply(1 - cp - mncp, 3:4, FUN = function(x) mymean(x > (1.2 * alpha / 2))), 3), nsmall = 1)
    pctCons.1side <- format(100 * round(apply(lncp, 3:4, FUN = function(x) mymean(x < (alpha / 2))), 3), nsmall = 1)
    pctnear.1side <- format(100 * round(apply(lncp, 3:4, FUN = function(x) mymean(x > 0.8 * (alpha / 2) & x < 1.2 * (alpha / 2))), 3), nsmall = 1)
    pctnear.DNCP <- format(100 * round(apply((1 - cp - mncp), 3:4, FUN = function(x) mymean(x > 0.8 * (alpha / 2) & x < 1.2 * (alpha / 2))), 3), nsmall = 1)
    pctCons.DNCP <- format(100 * round(apply((1 - cp - mncp), 3:4, FUN = function(x) mymean(x < (alpha / 2))), 3), nsmall = 1)
    pctBad.either <- format(100 * round(apply(maxncp, 3:4, FUN = function(x) mymean(x > (1.2 * alpha / 2))), 3), nsmall = 1)
    pctCons.both <- format(100 * round(apply(maxncp, 3:4, FUN = function(x) mymean(x < (alpha / 2))), 3), nsmall = 1)

    pctAveCons <-  pctAvenear <-  pctAvenearlo <- pctAvenearhi <- pctAveBad <-
      pctAvenear.1side <- typeI <- maxtypeI <- avetypeI <-
      array(NA, dim=c(max(length(psis), length(phis)), nmeth))

if (smooth == TRUE) {
#  option to add smoothed averages, but increases runtimes
  cat(paste0("Moving average coverage probabilities\n"))

      ## For a 'smoothed' surface plot, we can calculate the local average CP
      ## in the region (p1-avedelta, p1+avedelta, p2-avedelta, p2+avedelta)
      ma.win <- round(0.5 * window * length(p1))
      arraytemp <- array(NA, dim = dim(cp))
      dimnames(arraytemp) <- list(paste(p1), paste(p2), paste(par3), mymethods)
      avecp <- avelncp <- averncp <- avelocindex <- arraytemp

      # Calculate moving averages
#      pb <- txtProgressBar(min = 0, max = nmeth, style = 3) #text based bar
      pb <- pbapply::timerProgressBar(min = 0, max = nmeth, width = 50, char = "=", style = 1) #text based bar
      on.exit(pbapply::closepb(pb))
        for (meth in 1:nmeth) {
#          setTxtProgressBar(pb, meth)
          pbapply::setTimerProgressBar(pb, meth)
          for (p3 in par3) {
            avecp[-c(1:ma.win, (dim(cp)[1] - (ma.win - 1)):(dim(cp)[1])),
                  -c(1:ma.win, (dim(cp)[2] - (ma.win - 1)):(dim(cp)[2])),
                  paste(p3), meth] <-
              t(rollapply(zoo(t(rollapply(zoo(cp[, , paste(p3), meth]), (2 * ma.win + 1), mymean))), 2 * ma.win + 1, mymean))
            avelncp[-c(1:ma.win, (dim(cp)[1] - (ma.win - 1)):(dim(cp)[1])),
                    -c(1:ma.win, (dim(cp)[2] - (ma.win - 1)):(dim(cp)[2])),
                    paste(p3), meth] <-
              t(rollapply(zoo(t(rollapply(zoo(lncp[, , paste(p3), meth]), 2 * ma.win + 1, mymean))), 2 * ma.win + 1, mymean))
            averncp[-c(1:ma.win, (dim(cp)[1] - (ma.win - 1)):(dim(cp)[1])),
                    -c(1:ma.win, (dim(cp)[2] - (ma.win - 1)):(dim(cp)[2])),
                    paste(p3), meth] <-
              t(rollapply(zoo(t(rollapply(zoo(rncp[, , paste(p3), meth]), 2 * ma.win + 1, mymean))), 2 * ma.win + 1, mymean))
            avelocindex[-c(1:ma.win, (dim(cp)[1] - (ma.win - 1)):(dim(cp)[1])),
                        -c(1:ma.win, (dim(cp)[2] - (ma.win - 1)):(dim(cp)[2])),
                   paste(p3), meth] <-
              t(rollapply(zoo(t(rollapply(zoo(locindex[, , paste(p3), meth]), 2 * ma.win + 1, mymean))), 2 * ma.win + 1, mymean))
          }
        }
      pbapply::closepb(pb)

      arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "avecp", , ] <- avecp
      arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "avelncp", , ] <- avelncp
      arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "averncp", , ] <- averncp
      arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "avelocindex", , ] <- avelocindex
      # Note expected interval length tends to be smooth already

      pctAveCons <- format(100 * round(apply(avecp, 3:4, FUN = function(x) mean(x[!is.na(x)] > (1 - alpha), na.rm = T)), 3), nsmall = 1)
      pctAvenear <- format(100 * round(apply(avecp, 3:4, FUN = function(x) mean(x > (1 - 1.1 * alpha) & x < (1 - 0.9 * alpha), na.rm = T)), 3), nsmall = 1)
      pctAvenearlo <- format(100 * round(apply(avecp, 3:4, FUN = function(x) mean(x > (1 - 1.1 * alpha) & x < (1 - alpha), na.rm = T)), 3), nsmall = 1)
      pctAvenearhi <- format(100 * round(apply(avecp, 3:4, FUN = function(x) mean(x >= (1 - alpha) & x < (1 - 0.9 * alpha), na.rm = T)), 3), nsmall = 1)
      pctAveBad <- format(100 * round(apply(avecp, 3:4, FUN = function(x) mean(x <= (1 - alpha * 1.1), na.rm = T)), 3), nsmall = 1)
      pctAvenear.1side <- format(100 * round(apply(avelncp, 3:4, FUN = function(x) mean(x > 0.9 * (alpha / 2) & x < 1.1 * (alpha / 2), na.rm = T)), 3), nsmall = 1)

      rm(avecp, avelncp, averncp, avelocindex)

} # Close bracket for smoothing


    arrays[["summaries"]][,, paste(100 * (1 - alpha)), , 1, 1] <-
      array(c(meanCP, minCP, pctCons, pctBad, pctBad2, pctnear,
              pctAveCons, pctAvenear, pctAvenearlo, pctAvenearhi, pctAveBad,
              meanMNCP, meanDNCP, maxLNCP, pctCons.1side, pctBad.1side, pctBad.DNCP,
              pctnear.1side, pctAvenear.1side, typeI, maxtypeI, avetypeI,
              meanlen, meanlocindex, pctgoodloc,
              pctnear.DNCP, pctCons.DNCP, pctBad.either, pctCons.both),
            dim=c(max(length(psis), length(phis)), nmeth, length(summarylist)))

        # Free up some memory
    rm(cp, lncp, rncp, mncp, maxncp, locindex, len,
           cpl, lncpl, rncpl, mncpl, lenl, locindexl)

    save(arrays,
#         file = paste(outpath, "cparrays.", contrast, ".",
#                      n, ".", n.grid, ".Rdata", sep = "")
         file = paste("data/", "cparrays.", contrast, ".",
                      n, ".", n.grid, ".Rdata", sep = "")
    )
  }

  arrays
}


# Calculate coverage probabilities for a list (instead of a grid) of PSPs
# with pre-calculated ciarray, (for 2-d plots)
# or for a single PSP with calculation of CIs for all x's with non-negligible density
onecpfun <- function(
    p1 = c(0.3),
    p2 = c(0.1),
    ciarrays = NULL,
    n = NULL,
    contrast = NULL,
    alph = NULL,
    psis = NULL,
    phis = 0.25,
    methods = "All",
    jitt = FALSE
) {

  if (!is.null(ciarrays)) {
    xs <- ciarrays[["xs"]]
    cis <- ciarrays[["cis"]]
    contrast <- dimnames(cis)[[6]]
    mymethods <- longlab <- dimnames(cis)[[3]]
    nmeth <- length(mymethods)
    n <- as.numeric(dimnames(cis)[[5]])

  if (!is.null(psis)) par3 <- psis
  if (!is.null(phis)) par3 <- rep_len(phis, length.out = length(p1))

  mastercp <- array(NA, dim = c(length(p1), length(mymethods), length(alph), 6, 1, 1))
  dimnames(mastercp) <- list(paste(p2),
                             mymethods,
                             paste(100*(1-alph)),
                             c("cp", "lncp", "rncp", "dncp", "len", "locindex"),
                             paste(n),
                             contrast)

  px <- pxg <- cbind(p1, p2, par3)

  if (contrast %in% c("RR")) {
    theta <- (px[, 1] / px[, 2])
  } else if (contrast == "OR") {
    theta <- (px[, 1] * (1 - px[, 2]) / (px[, 2] * (1 - px[, 1])))
  } else if (contrast == "RD") {
    theta <- px[, 1] - px[, 2]
  } # theta estimate for every p1,p2 pair

  cat(paste0("\nStep 2 of 3: calculate coverage probabilities for each alpha"))

  for (alpha in alph) {
    ci <- cis[,,, paste(100 * (1 - alpha)),,]

    # To obtain coverage probability for n,p1,p2, we determine the probability
    # of each possible pair of observed frequencies a,b given p1,p2 using
    # binomial probabilities
    # and sum over all combinations where p1-p2 is included in the corresponding
    #  CI, for each CI method.
    cpl <- lncpl <- rncpl <- lenl <- locindexl <-
      array(NA, dim = c(dim(px)[1], nmeth))

    # NB method using sapply was no quicker than this loop:
    # NB this is the most time-consuming part of the process when n is small
    # i <- 40 #gives negative probs at p1=0.95, p2=0.35. phi=0.5
    # - pdfpair function updated to give prob=NA for such impossible combinations
    # First switched to parameterisation using psi instead.
    # Then revisited in order to produce plots based on phi (more interpretable)
    cat(paste0("\nalpha=", alpha))
#    pb <- txtProgressBar(min = 0, max = dim(px)[1], style = 3) #text based bar
    pb <- pbapply::timerProgressBar(min = 0, max = dim(px)[1], width = 50, char = '=', style = 1) #text based bar
    on.exit(pbapply::closepb(pb))
    for (i in 1:dim(px)[1]) {
#      setTxtProgressBar(pb, i)
      pbapply::setTimerProgressBar(pb, i)
      if (!is.null(psis)) {
        prob <- pdfpair(p1 = px[i, 1],
                        p2 = px[i, 2],
                        psi = px[i, 3],
                        x = xs)
      } else if (!is.null(phis)) {
        prob <- pdfpair(p1 = px[i, 1],
                        p2 = px[i, 2],
                        phi = px[i, 3],
                        x = xs)
      }

      if (any(!is.na(prob))) {

        cpl[i, ] <- t(ci[, 1, ] <= theta[i] & ci[, 2, ] >= theta[i] & ci[, 2, ] > ci[, 1, ]) %*% prob # 2-sided coverage probability. NB degenerate intervals excluded
        lncpl[i, ] <- t(ci[, 1, ] > theta[i] | ci[, 2, ] == ci[, 1, ]) %*% prob # L-sided non-coverage (R-side is a mirror image)
        rncpl[i, ] <- t(ci[, 2, ] < theta[i] | ci[, 2, ] == ci[, 1, ]) %*% prob # R-sided non-coverage

        # not really interested in length, I think location is more important
        # - and its complicated to explain on log scale
        # but included as reviewers might request it
        lens <- ci[, 2, ] - ci[, 1, ]
        ## you get  whacky results with length on linear scale
        # 		if(contrast %in% c("RR")) lens<-(ci[,2,])-(ci[,1,])
        # 		lens[lens>100]<-100 #workaround for infinite lengths with RR
        if (contrast %in% c("RR")) {
          # for RR, use length on the log scale
      #    lens <- log(ci[, 2, ]) - log(pmax(0.0000000001, ci[, 1, ]))
      #    lens <- ci[, 2, ]/(1 + ci[, 2, ]) - ci[, 1, ]/(1 + ci[, 1, ])
          lens <- log(1/(1 + 1/ci[, 2, ])) - log(ci[, 1, ]/(1 + ci[, 1, ]))
        }
#        lens[lens > 10] <- 10 # workaround for infinite lengths with RR
        lens[ci[, 2, ] == ci[, 1, ]] <- 0
        lenl[i, ] <- t(lens) %*% prob
      }
    }

    mncpl <- rncpl
    thetagt <- px[,1] > px[,2]
    mncpl[thetagt, ] <- lncpl[thetagt, ]
    locindexl <- ifelse(cpl == 1, 0, mncpl / (1 - cpl))
    dncpl <- 1 - cpl - mncpl

    # Free up some memory
    rm(ci)

    dimnames(cpl) <- dimnames(locindexl) <- dimnames(lncpl) <-
      dimnames(rncpl) <- dimnames(lenl) <- dimnames(dncpl) <-
      list(p2, mymethods)

    mastercp <- array(rbind(cpl, lncpl, rncpl, dncpl, locindexl, lenl), dim = c(length(p2), 6, nmeth))
    dimnames(mastercp) <- list(p2, c("cp", "lncp", "rncp", "dncp", "locindex", "len") , mymethods)

    # Free up some memory
    rm(
       cpl, lncpl, rncpl, dncpl, lenl, locindexl)

  }

  myoutput <- aperm(mastercp, c(1,3,2))

  save(mastercp,
       file = paste(outpath, "mastercp.", contrast, ".",
                    n, ".Rdata", sep = "")
  )

}

  ### Following section runs calculations for a single PSP ###
  ### without having pre-run the CIs ###
  if (!is.null(n) && length(p1) == 1 && length(p2) == 1 && length(phis) == 1) {
    nmeth <- length(methods)
    # For a given N, find the set of outcomes with non-negligible probabilities
    g <- (expand.grid(x11 = 0:n, x12 = 0:n, x21 = 0:n))
    # reduce to possible combinations of a,b,c for paired data.
    g <- g[(g$x12 <= n - g$x11) &
             (g$x21 <= n - g$x11 - g$x12), ]
    xs <- data.matrix(cbind(g, x22 = n - rowSums(g)))
    lenxs <- dim(xs)[1]
    row.names(xs) <- NULL
    px <- array(c(p1, p2, phis), dim = c(1, 3))
    i <- 1 # For now, only allow this for a single PSP
    if (!is.null(psis)) {
      prob <- pdfpair(p1 = px[i, 1],
                      p2 = px[i, 2],
                      psi = px[i, 3],
                      x = xs)
    } else if (!is.null(phis)) {
      prob <- pdfpair(p1 = px[i, 1],
                      p2 = px[i, 2],
                      phi = px[i, 3],
                      x = xs)
    }

    xsub <- xs[prob > 1E-10, ]
    cpl <- lncpl <- rncpl <- lenl <- locindexl <-
      array(NA, dim = c(dim(px)[1], nmeth))

    tester <- allpairci(x = rep(10,4), contrast = contrast)
    if (is.null(methods)) {
      mymethods <- dimnames(tester)[[3]]
    } else mymethods <- methods
    nmeth <- length(mymethods)

#    ci <- allpairci(x = xsub, contrast = contrast, alpha = alph, methods = mymethods)
    ci <- array(NA, dim=c(dim(xsub)[1], 2, nmeth))
    # Process one x at a time so we can show progress and expected time
    pb <- pbapply::timerProgressBar(min = 0, max = dim(xsub)[1], width = 50, char = '=', style = 1) #timer bar
    on.exit(pbapply::closepb(pb))
    for (j in 1:dim(xsub)[[1]]) {
      pbapply::setTimerProgressBar(pb, j)
      ci1 <- allpairci(xs = xsub[j, ], contrast = contrast, alpha = alph, methods = mymethods)
      ci[j, , ] <- ci1
    }
    pbapply::closepb(pb)

    if (contrast == "RD") theta <- p1 - p2
    if (contrast == "RR") theta <- p1 / p2
    if (contrast == "OR") {
      pij <- params(p1 = p1, p2 = p2, phi = phis, psi = psis)[, 5:8]
      theta <- pij[2] / pij[3] # Conditional OR
    }

    cpl[i, ] <- t(ci[, 1, ] <= theta[i] & ci[, 2, ] >= theta[i] & ci[, 2, ] > ci[, 1, ]) %*% prob[prob > 1E-10] # 2-sided coverage probability. NB degenerate intervals excluded
    lncpl[i, ] <- t(ci[, 1, ] > theta[i] | ci[, 2, ] == ci[, 1, ]) %*% prob[prob > 1E-10] # L-sided non-coverage (R-side is a mirror image)
    rncpl[i, ] <- t(ci[, 2, ] < theta[i] | ci[, 2, ] == ci[, 1, ]) %*% prob[prob > 1E-10] # R-sided non-coverage

    mncpl <- rncpl
    thetagt <- px[,1] > px[,2]
    mncpl[thetagt, ] <- lncpl[thetagt, ]
    locindexl <- ifelse(cpl == 1, 0, mncpl / (1 - cpl))
    dncpl <- 1 - cpl - mncpl

    mastercp <- array(rbind(cpl, lncpl, rncpl, dncpl, locindexl, lenl), dim = c(length(p2), 6, nmeth))
    dimnames(mastercp) <- list(p2, c("cp", "lncp", "rncp", "dncp", "locindex", "len") , mymethods)

    myoutput <- aperm(mastercp, c(3,2,1))

  }

  return(myoutput)
}

