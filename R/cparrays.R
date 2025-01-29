# Evaluation of confidence intervals for RD or RR from paired data
# 2025 update using parameter p11p22/p12p21 from Fagerland instead of rho(=phi)
# and adding MOVER methods
#install.packages("zoo")

rm(list = ls())
library(zoo)
pak::pak("petelaud/ratesci-dev")
library(ratesci)

set.seed(2012) #ensure we use the same jitters for each run

root <- "/Users/ssu/Documents/"
outpath <- paste(root, "Main/Courses_papers/skewscore/paired/", sep = "")

#source(paste(newpath,'diffbinplotsv4.paired.R',sep=""))

mysum <- function(x) sum(x,na.rm=T)
mymean <- function(x) mean(x,na.rm=T)

# Conversion of parameters for RR
params <- function(p1,
                   theta,
                   psi=NULL,
                   phi=NULL
                   ) {
  # convert phi to psi, or vice versa
  p2 <- p1 / theta
  if (is.null(phi)) {
    A <- psi-1
    B <- p1+p2-1-psi*(p1+p2)
    C <- psi*p1*p2
    p11 <- (-B - sqrt(B^2-4*A*C))/(2*A)
  }
  if (is.null(psi)) {
    p11 <- p1 * p2 + phi*sqrt(p1*(1-p1)*p2*(1-p2)) # with DelRocco parameterisation using rho(=phi)
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
par(mfrow=c(1,3), pty='s', mar=(c(2,4,2,0)+0.1))
plot(p1s, params(p1=p1s, theta=4, psi=2)[,4], ylim = c(0, 0.9), type='n',
     main=expression(paste(theta["RR"], " = 1")), xlab = "p1",
     ylab = expression(phi))
lines(p1s, params(p1=p1s, theta=1, psi=2)[,4], lty=1)
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

# For a given parameter combination: p1, p2, (psi or phi)
# calculate the probability density function for each combination
pdfpair <- function(x,
                    p1 = NULL,
                    p2 = NULL,
                    phi = NULL,
                    psi = NULL
                    ) {
  if (is.null(psi)) {
    p11 <- p1 * p2 + phi * sqrt(p1 * (1 - p1) * p2 * (1 - p2)) # with DelRocco parameterisation using rho(=phi)
  }
  # Or Fagerland parameterisation using psi=(p11*p22)/(p12*p21) (labelled as theta in Fagerland)
  if (is.null(phi)) {
    A <- psi - 1
    B <- p1 + p2 - 1 - psi * (p1 + p2)
    C <- psi * p1 * p2
    p11 <- ifelse(psi == 1, -C / B,(-B - sqrt(B^2 - 4 * A * C)) / (2 * A))
  }
  p12 <- pmax(0, p1 - p11)
  p21 <- pmax(0, p2 - p11)
  p22 <- pmax(0, 1 - p11 - p12 - p21)
  dens <- exp(
    lfactorial(rowSums(x)) -
      (lfactorial(x[,1]) + lfactorial(x[,2]) + lfactorial(x[,3]) + lfactorial(x[,4]))
  ) * p11^x[,1] * p12^x[,2] * p21^x[,3] * (1 - p11 - p12 - p21)^(x[,4])
  dens[p21 == 0 | p12 == 0] <- NA
  return(dens)
}

cifun <- function( n,
                   contrast = "RD",
                   alph = c(0.01, 0.05, 0.1),
                   prerun = FALSE
                   ) {
  # This function calculates all possible confidence intervals for a given n,

  if (FALSE) {
  # test run code
  # m=15; n=15; alph=0.05; p2=0.4; p1=0.4; exact=F; contrast="RR"; window=0.08; square=T; prerun=F; jitt=F
  # alpha=0.05; sided="R"
  # n=15
  # i=23
  # 	limits <- c(0,1); n.grid <- 10
  # 	refine=1/n.grid
  # 	p.grid<-limits[1]+ppoints(n.grid,a=0.5)*(limits[2]-limits[1])
  # 	p1<-p2<-p.grid #+jitter #add random jitter here if required
  # 	phi=0

#  n <- 20
#  contrast <- "RD"
  }

  if (contrast == "RR") {
    mymethods <- c("Tang", "Tang-bc", "Tang-sc", "Tang-scbc",
                   "MOVER-w", "MOVER-j", "MOVER-s",
                   "MOVER-nw", "MOVER-nj", "MOVER-ns", "Tang-cc5",
                   "Tang-ccdr", "Tang-cc25", "Tang-cc125", "Tang-sccc",
                   "MOVER-ccnw", "MOVER-ccns", "TDAS")
  } else if (contrast == "RD") {
    mymethods <- c("Tango", "Tango-bc", "Tango-sc", "Tango-scbc",
                   "MOVER-w", "MOVER-j", "MOVER-s",
                   "MOVER-nw", "MOVER-nj", "MOVER-ns", "Tango-cc5",
                   "Tango-cc25", "Tango-cc125", "Tango-sccc",
                   "MOVER-ccnw", "MOVER-ccns", "TDAS")
  } else if (contrast == "OR") {
    mymethods <- c("SCAS", "Jeffreys", "mid-p", "Wilson")
  }

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

  cis <- array(NA, dim = c(lenxs, 2, length(mymethods), length(alph), 1, 1))
  dimnames(cis)[[2]] <- c("LCL", "UCL")
  dimnames(cis)[[3]] <- mymethods
  dimnames(cis)[[4]] <- paste(100*(1-alph))
  dimnames(cis)[[5]] <- paste(n)
  dimnames(cis)[[6]] <- contrast

  for (alpha in alph) {
    # For the given sample size m, we calculate confidence intervals,
    # for each method, for every possible combination of observed frequencies

    # Might be able to reuse same code for RD & RR?
    if (contrast == "RR") {
        ci <- array(NA, dim = c(lenxs, 2, length(mymethods)))
        dimnames(ci)[[3]] <- mymethods
        ci[, 1:2, "Tang"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "Tang-bc"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RR = "Score_closed", bcf=TRUE, level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "Tang-sc"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RR = "Score", skew=TRUE, level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "Tang-scbc"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RR = "Score", skew=TRUE, bcf=TRUE, level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "Tang-ccdr"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, level = 1-alpha, cc=0.5, cctype="delrocco")$estimates[,c(1,3)]))
        ci[, 1:2, "Tang-cc5"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, level = 1-alpha, cc=0.5, cctype="new")$estimates[,c(1,3)]))
        ci[, 1:2, "Tang-cc25"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, level = 1-alpha, cc=0.25, cctype="new")$estimates[,c(1,3)]))
        ci[, 1:2, "Tang-cc125"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, level = 1-alpha, cc=0.125, cctype="new")$estimates[,c(1,3)]))
        ci[, 1:2, "Tang-sccc"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RR = "Score", skew=TRUE, cc=0.125, cctype="new", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-w"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RR = "MOVER", moverbase = "wilson", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-j"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RR = "MOVER", moverbase = "jeff", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-s"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RR = "MOVER", moverbase = "SCAS", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-nw"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RR = "MOVER_newc", moverbase = "wilson", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-nj"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RR = "MOVER_newc", moverbase = "jeff", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-ns"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RR = "MOVER_newc", moverbase = "SCAS", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-ccnw"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RR = "MOVER_newc", moverbase = "wilson", level = 1-alpha, cc=0.125)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-ccns"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RR = "MOVER_newc", moverbase = "SCAS", level = 1-alpha, cc=0.125)$estimates[,c(1,3)]))
        # TDAS takes a while - leave out for now
#        if (alpha == 0.05) {
#          ci[, 1:2, "TDAS"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RR = "TDAS")$estimates[c(1,3)]))
#        }
    } else if (contrast == "RD") {
        ci <- array(NA, dim = c(lenxs, 2, length(mymethods)))
        dimnames(ci)[[3]] <- mymethods
        ci[, 1:2, "Tango"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "Tango-bc"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, bcf = TRUE, level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "Tango-sc"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RD = "Score", skew=TRUE, level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "Tango-scbc"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RD = "Score", skew=TRUE, bcf=TRUE, level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "Tango-cc125"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, cc = 0.125, level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "Tango-cc25"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, cc = 0.25, level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "Tango-cc5"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, cc = 0.5, level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "Tango-sccc"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RD = "Score", skew=TRUE, cc=0.125, level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-w"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RD = "MOVER", moverbase = "wilson", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-j"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RD = "MOVER", moverbase = "jeff", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-s"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RD = "MOVER", moverbase = "SCAS", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-nw"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RD = "MOVER_newc", moverbase = "wilson", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-nj"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RD = "MOVER_newc", moverbase = "jeff", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-ns"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RD = "MOVER_newc", moverbase = "SCAS", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-ccnw"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RD = "MOVER_newc", moverbase = "wilson", cc=0.125, level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "MOVER-ccns"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RD = "MOVER_newc", moverbase = "SCAS", cc=0.125, level = 1-alpha)$estimates[,c(1,3)]))
        # TDAS takes a while - leave out for now
#        if (alpha == 0.05) {
#          ci[, 1:2, "TDAS"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_RD = "TDAS")$estimates[c(1,3)]))
#        }
    } else if (contrast == "OR") {
        ci <- array(NA, dim = c(lenxs, 2, length(mymethods)))
        dimnames(ci)[[3]] <- mymethods
        ci[, 1:2, "SCAS"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_OR = "SCAS", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "Jeffreys"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_OR = "jeff", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "mid-p"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_OR = "midp", level = 1-alpha)$estimates[,c(1,3)]))
        ci[, 1:2, "Wilson"] <- t(sapply(1:lenxs,function(i) pairbinci(x = xs[i,], contrast = contrast, method_OR = "wilson", level = 1-alpha)$estimates[,c(1,3)]))
    }

    cis[, , , paste(100 * (1 - alpha)), 1, 1] <- ci

  }

if (FALSE) {
  p.grid <- limits[1] + ppoints(n.grid, a = 0.5) * (limits[2] - limits[1])
  nn <- qbinom(0.99999, n, limits[2]) # Upper limit for n when considering Poisson rates?
  p1 <- p2 <- p.grid

  if (!is.null(psis)) par3 <- psis
  if (!is.null(phis)) par3 <- phis

  mastercp <- array(NA, dim = c(length(p1), length(p2), length(par3), length(mymethods), length(alph), 10, 1))
  dimnames(mastercp) <- list(paste(p1),
                             paste(p2),
                             paste(par3),
                             mymethods,
                             paste(100*(1-alph)),
                             c("cp", "lncp", "rncp", "mncp", "avecp",
                               "avelncp", "averncp", "len", "locindex", "avelocindex"),
                             paste(n))

  summarylist <- c("meanCP", "minCP", "pctCons", "pctBad", "pctnear",
    "pctAveCons", "pctAvenear", "pctAvenearlo", "pctAvenearhi", "pctAveBad",
    "meanMNCP", "meanDNCP", "maxLNCP", "pctCons.1side", "pctBad.1side",
    "pctnear.1side", "pctAvenear.1side", "typeI", "maxtypeI", "avetypeI",
    "meanlen", "meanlocindex", "pctgoodloc")
  summaries <- array(NA, dim = c(length(par3), length(mymethods), length(alph), length(summarylist), 1))
  dimnames(summaries) <-
    list(paste(par3),
         mymethods,
         paste(100*(1-alph)),
         summarylist,
         paste(n))

#  if (prerun == FALSE) {
  myarrays <- list(xs = xs, cis = cis, mastercp = mastercp, summaries = summaries)
  rm(mastercp)
  rm(summaries)
  rm(cis)

  save(myarrays,
       file = paste(outpath, "arrays.",
                    ifelse(contrast == "OR", "OR.",
                           ifelse(contrast == "RR", "RR.",
                                  ifelse(contrast == "ID", "ID.",
                                         ifelse(contrast == "IDRR", "IDRR.", "")))),
                    n, ".", n.grid, ".Rdata", sep = "")
  )
}
  ciarrays <- list(xs = xs, cis = cis)
  save(ciarrays, file = paste0(outpath, "cis.", contrast, ".", n, ".Rdata"))
  ciarrays
}

#myciarrays <- cifun(n=10, contrast="RD", alph=c(0.05, 0.1))

#    arrays = myarrays,
cpfun <- function(
                  ciarrays = myciarrays,
#                  contrast = "RR",
#                  alph = NULL,
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
                  sided = "R"
                  ) {
  # This function calculates coverage probabilities (raw and smoothed) and associated summaries

  xs <- ciarrays[["xs"]]
  cis <- ciarrays[["cis"]]
  contrast <- dimnames(cis)[[6]]
  mymethods <- longlab <- dimnames(cis)[[3]]
  alph <- 1 - (as.numeric(dimnames(cis)[[4]])/100)
  nmeth <- length(mymethods)
  n <- as.numeric(dimnames(cis)[[5]])

#  p1 <- as.numeric(dimnames(arrays[["mastercp"]])[[1]])
#  p2 <- as.numeric(dimnames(arrays[["mastercp"]])[[2]])
#  if (is.null(psis)) psis <- as.numeric(dimnames(arrays[["mastercp"]])[[3]])
#  if (is.null(alph)) alph <- (1 - as.numeric(dimnames(arrays[["mastercp"]])[[5]])/100)
#    names(longlab) <- dimnames(cis)[[3]]

  if (!is.null(psis)) par3 <- psis
  if (!is.null(phis)) par3 <- phis
  p.grid <- limits[1] + ppoints(n.grid, a = 0.5) * (limits[2] - limits[1])
  nn <- qbinom(0.99999, n, limits[2]) # Upper limit for n when considering Poisson rates?
  p1 <- p2 <- p.grid

  mastercp <- array(NA, dim = c(length(p1), length(p2), length(par3), length(mymethods), length(alph), 10, 1, 1))
  dimnames(mastercp) <- list(paste(p1),
                             paste(p2),
                             paste(par3),
                             mymethods,
                             paste(100*(1-alph)),
                             c("cp", "lncp", "rncp", "mncp", "avecp",
                               "avelncp", "averncp", "len", "locindex", "avelocindex"),
                             paste(n),
                             contrast)

  summarylist <- c("meanCP", "minCP", "pctCons", "pctBad", "pctnear",
                   "pctAveCons", "pctAvenear", "pctAvenearlo", "pctAvenearhi", "pctAveBad",
                   "meanMNCP", "meanDNCP", "maxLNCP", "pctCons.1side", "pctBad.1side",
                   "pctnear.1side", "pctAvenear.1side", "typeI", "maxtypeI", "avetypeI",
                   "meanlen", "meanlocindex", "pctgoodloc")
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
  set.seed(2012) # ensure we use the same jitters for each run
  px <- pxg <- expand.grid(p1, p2, par3)

  refine <- 1 / n.grid
  if (jitt) {
    jitter <- (limits[2] - limits[1]) * array(runif(2 * dim(pxg)[1], -refine / 2, refine / 2), dim = dim(pxg)) # OK if n.grid is even
    px[, 1:2] <- pxg[, 1:2] + jitter[, 1:2]
  } else {
    px <- pxg
  }

  # dimnames(px)[2]=c("p1","p2")
  if (contrast %in% c("RR")) {
    theta <- (px[, 1] / px[, 2])
  } else if (contrast == "OR") {
#    theta <- log(px[, 1] * (1 - px[, 2]) / (px[, 2] * (1 - px[, 1])))
    theta <- (px[, 1] * (1 - px[, 2]) / (px[, 2] * (1 - px[, 1])))
  } else if (contrast == "RD") {
    theta <- px[, 1] - px[, 2]
  } # theta estimate for every p1,p2 pair

  # can we improve efficiency here? -
  # not without creating an array which is too large for available RAM.

# alpha <- 0.05
  for (alpha in alph) {
#    if (contrast != "OR") {
      ci <- cis[,,, paste(100 * (1 - alpha)),,]
#    } else ci <- log(cis[,,, paste(100 * (1 - alpha)),])

    # To obtain coverage probability for n,p1,p2, we determine the probability
    # of each possible pair of observed frequencies a,b given p1,p2 using
    # binomial probabilities
    # and sum over all combinations where p1-p2 is included in the corresponding
    #  CI, for each CI method.
    cpl <- lncpl <- rncpl <- lenl <- locindexl <-
      array(NA, dim = c(dim(px)[1], nmeth))

    # NB method using sapply was no quicker than this loop:
    # NB this is the most time-consuming part of the process when n is small
    # i <- 40 #gives negative probs at p1=0.95, p2=0.35. phi=0.5 - need to produce a warning for impossible parameter combinations
    # so I've switched to parameterisation using psi instead.
    for (i in 1:dim(px)[1]) {
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
      # - could loop by px and then by alpha?

      # exclude (0,0) CI from CP calcs by setting its probability to zero and
      # rescaling all other probabilities - not advisable
      # prob<-prob*(1-prob[1])
      # prob[1]<-0

#      cbind(xs, prob[, i])
#      prob[prob == 0] <- NA

      cpl[i, ] <- t(ci[, 1, ] <= theta[i] & ci[, 2, ] >= theta[i] & ci[, 2, ] > ci[, 1, ]) %*% prob # 2-sided coverage probability. NB degenerate intervals excluded
#        cpl[i, ] <-  t(prob) %*% (ci[, 1, ] <= theta[i] & ci[, 2, ] >= theta[i] & ci[, 2, ] > ci[, 1, ])
      lncpl[i, ] <- t(ci[, 1, ] > theta[i] | ci[, 2, ] == ci[, 1, ]) %*% prob # L-sided non-coverage (R-side is a mirror image)
      rncpl[i, ] <- t(ci[, 2, ] < theta[i] | ci[, 2, ] == ci[, 1, ]) %*% prob # R-sided non-coverage

if (FALSE) {
# Exploration of OR issues
#      rmiss <- ci[, 2, ] < theta[i] | ci[, 2, ] == ci[, 1, ]
#      cbind(xs,rmiss,ci[,2,], prob)[rowSums(rmiss) > 0 & rowSums(rmiss) < 4 & prob > 1e-2, ]
#       Jeffreys misses c(6,0,7,7)
#       True OR = 0.444 p1*(1-p2)/(p2*(1-p1))
#      pairbinci(x=c(6,0,7,7), contrast="OR", method_OR="jeff")
#      pairbinci(x=c(6,0,7,7), contrast="OR", method_OR="SCAS")
#      pairbinci(x=c(6,0,7,7), contrast="OR", method_OR="wilson")
# Margin where rncp=0 is bigger because the conditional method is working with less data?
}

      # not really interested in length, I think location is more important
      # - and its complicated to explain on log scale
      # but included as reviewers might request it
      lens <- ci[, 2, ] - ci[, 1, ]
      ## you get pretty whacky results with length on linear scale
      # 		if(contrast %in% c("RR")) lens<-(ci[,2,])-(ci[,1,])
      # 		lens[lens>100]<-100 #workaround for infinite lengths with RR
      if (contrast %in% c("RR")) {
        # for RR, use length on the log scale
        lens <- log(ci[, 2, ]) - log(pmax(0.0000000001, ci[, 1, ]))
      }
      lens[lens > 10] <- 10 # workaround for infinite lengths with RR
      lens[ci[, 2, ] == ci[, 1, ]] <- 0
      lenl[i, ] <- t(lens) %*% prob
    }

    mncpl <- rncpl
    thetagt <- px[,1] > px[,2]
    mncpl[thetagt, ] <- lncpl[thetagt, ]

    # Free up some memory
    rm(ci)

    # convert CPs etc to a square grid of p1,p2 points
    mydim <- c(length(p1), length(p2), length(par3), nmeth)
    cp <- array(cpl, dim = mydim)
    lncp <- array(lncpl, dim = mydim)
    rncp <- array(rncpl, dim = mydim)
    mncp <- array(mncpl, dim = mydim)
    len <- array(lenl, dim = mydim)
    # Calculate Newcombe's location index
    locindex <- ifelse(cp == 1, 0, mncp / (1 - cp))
    dimnames(cp) <- dimnames(locindex) <- dimnames(lncp) <-
      dimnames(rncp) <- dimnames(len) <- dimnames(mncp) <-
      list(paste(p1), paste(p2), paste(par3), mymethods)

    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "cp", , ] <- cp
    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "lncp", , ] <- lncp
    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "rncp", , ] <- rncp
#    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "mncp", , ] <- mncp
    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "locindex", , ] <- locindex
    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "len", , ] <- len

    meanlen <- format(round(apply(len, 3:4, mymean), 3), nsmall = 3)
    meanCP <- format(round(apply(cp, 3:4, mymean), 3), nsmall = 3)
    minCP <- format(round(apply(cp, 3:4, function(x) min(x, na.rm = TRUE)), 3), nsmall = 3)
    pctCons <- format(100 * round(apply(cp, 3:4, FUN = function(x) mymean(x > (1 - alpha))), 3), nsmall = 1)
    pctBad <- format(100 * round(apply(cp, 3:4, FUN = function(x) mymean(x < (1 - (alpha * 1.1)))), 3), nsmall = 1)
    pctnear <- format(100 * round(apply(cp, 3:4, FUN = function(x) mymean(x > (1 - 1.1 * alpha) & x < (1 - 0.9 * alpha))), 3), nsmall = 1)
    meanlocindex <- format(round(apply(locindex, 3:4, mymean), 3), nsmall = 3)
    pctgoodloc <- format(100 * round(apply(locindex, 3:4, FUN = function(x) mymean(x > 0.4 & x < 0.6)), 3), nsmall = 1)
    # one-sided summaries
    meanMNCP <- format(round(apply(mncp, 3:4, mymean), 3), nsmall = 3)
    meanDNCP <- format(round(apply((1 - cp - mncp), 3:4, mymean), 3), nsmall = 3)
    maxLNCP <- format(round(apply(lncp, 3:4, function(x) max(x, na.rm = TRUE)), 3), nsmall = 3)
    pctBad.1side <- format(100 * round(apply(lncp, 3:4, FUN = function(x) mymean(x > (1.1 * alpha / 2))), 3), nsmall = 1)
    pctCons.1side <- format(100 * round(apply(lncp, 3:4, FUN = function(x) mymean(x < (alpha / 2))), 3), nsmall = 1)
    pctnear.1side <- format(100 * round(apply(lncp, 3:4, FUN = function(x) mymean(x > 0.9 * (alpha / 2) & x < 1.1 * (alpha / 2))), 3), nsmall = 1)

    pctAveCons <-  pctAvenear <-  pctAvenearlo <- pctAvenearhi <- pctAveBad <-
      pctAvenear.1side <- typeI <- maxtypeI <- avetypeI <-
      array(NA, dim=c(max(length(psis), length(phis)), nmeth))

if (smooth == TRUE) {
#  option to add smoothed averages, but increases runtimes substantially

      ## For a 'smoothed' surface plot, we can calculate the local average CP
      ## in the region (p1-avedelta, p1+avedelta, p2-avedelta, p2+avedelta)
      ma.win <- round(0.5 * window * length(p1))
      arraytemp <- array(NA, dim = dim(cp))
      dimnames(arraytemp) <- list(paste(p1), paste(p2), paste(par3), mymethods)
      avecp <- avelncp <- averncp <- avelocindex <- arraytemp
      if (square == TRUE) { # skip this bit if we're not wanting the whole square grid (usually do)
        pb <- txtProgressBar(min = 0, max = nmeth, style = 3) #text based bar
        for (meth in 1:nmeth) {
          setTxtProgressBar(pb, meth)
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
      }

      arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "avecp", , ] <- avecp
      arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "avelncp", , ] <- avelncp
      arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "averncp", , ] <- averncp
      arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "avelocindex", , ] <- avelocindex
      # Note expected interval length tends to be smooth already

if (FALSE) {
  # For reference
      c("meanCP", "minCP", "pctCons", "pctBad", "pctBad2", "pctnear",
        "pctAveCons", "pctAvenear", "pctAvenearlo", "pctAvenearhi", "pctAveBad",
        "meanMNCP", "meanDNCP", "maxLNCP", "pctCons.1side", "pctBad.1side",
        "pctnear.1side", "pctAvenear.1side", "typeI", "maxtypeI", "avetypeI",
        "meanlen", "meanlocindex")
}

      pctAveCons <- format(100 * round(apply(avecp, 3:4, FUN = function(x) mean(x[!is.na(x)] > (1 - alpha), na.rm = T)), 3), nsmall = 1)
      pctAvenear <- format(100 * round(apply(avecp, 3:4, FUN = function(x) mean(x > (1 - 1.1 * alpha) & x < (1 - 0.9 * alpha), na.rm = T)), 3), nsmall = 1)
      pctAvenearlo <- format(100 * round(apply(avecp, 3:4, FUN = function(x) mean(x > (1 - 1.1 * alpha) & x < (1 - alpha), na.rm = T)), 3), nsmall = 1)
      pctAvenearhi <- format(100 * round(apply(avecp, 3:4, FUN = function(x) mean(x >= (1 - alpha) & x < (1 - 0.9 * alpha), na.rm = T)), 3), nsmall = 1)
      pctAveBad <- format(100 * round(apply(avecp, 3:4, FUN = function(x) mean(x <= (1 - alpha * 1.1), na.rm = T)), 3), nsmall = 1)
      pctAvenear.1side <- format(100 * round(apply(avelncp, 3:4, FUN = function(x) mean(x > 0.9 * (alpha / 2) & x < 1.1 * (alpha / 2), na.rm = T)), 3), nsmall = 1)
#      typeI <- format(round(apply(lncp, 3:4, function(x) mean(diag(x[(n.grid / 2 + 1):(n.grid * 0.85), (n.grid * 0.6 + 1):(n.grid * 0.95)]))), 4), nsmall = 4) # Average LNCP along the delta=-0.1 diagonal, from p1=0.5 to 0.9
#      maxtypeI <- format(round(apply(lncp, 3:4, function(x) max(diag(x[(n.grid / 2 + 1):(n.grid * 0.85), (n.grid * 0.6 + 1):(n.grid * 0.95)]))), 4), nsmall = 4) # Average LNCP along the delta=-0.1 diagonal, from p1=0.5 to 0.9
#      avetypeI <- format(round(apply(avelncp, 3:4, function(x) mean(diag(x[(n.grid / 2 + 1):(n.grid * 0.85), (n.grid * 0.6 + 1):(n.grid * 0.95)]), na.rm = T)), 4), nsmall = 4) # Average LNCP along the delta=-0.1 diagonal, from p1=0.5 to 0.9

      rm(avecp, avelncp, averncp, avelocindex)

# Close bracket for smoothing
}

    arrays[["summaries"]][,, paste(100 * (1 - alpha)), , 1, 1] <-
      array(c(meanCP, minCP, pctCons, pctBad, pctBad2, pctnear,
              pctAveCons, pctAvenear, pctAvenearlo, pctAvenearhi, pctAveBad,
              meanMNCP, meanDNCP, maxLNCP, pctCons.1side, pctBad.1side,
              pctnear.1side, pctAvenear.1side, typeI, maxtypeI, avetypeI,
              meanlen, meanlocindex, pctgoodloc),
            dim=c(max(length(psis), length(phis)), nmeth, length(summarylist)))

        # Free up some memory
    rm(cp, lncp, rncp, mncp, locindex, len,
           cpl, lncpl, rncpl, mncpl, lenl, locindexl)

    save(arrays,
         file = paste(outpath, "cparrays.", contrast, ".",
                      n, ".", n.grid, ".Rdata", sep = "")
    )
  }

  arrays
}



onecpfun <- function(
    p1 = c(0.3, 0.4),
    p2 = c(0.1, 0.2),
    ciarrays = myciarrays,
    alph = NULL,
    psis = NULL,
    phis = 0.25,
#    n.grid = 200,
#    res.factor = 6,
#    methods = "All",
#    limits = c(0, 1),
#    window = 0.08,
#    square = TRUE,
#    smooth = TRUE,
#    prerun = F,
    jitt = FALSE
#    sided = "R"
) {
  # This function calculates coverage probabilities (raw and smoothed) and associated summaries

  xs <- ciarrays[["xs"]]
  cis <- ciarrays[["cis"]]
  contrast <- dimnames(cis)[[6]]
  mymethods <- longlab <- dimnames(cis)[[3]]
#  alph <- 1 - (as.numeric(dimnames(cis)[[4]])/100)
  nmeth <- length(mymethods)
  n <- as.numeric(dimnames(cis)[[5]])

  #  p1 <- as.numeric(dimnames(arrays[["mastercp"]])[[1]])
  #  p2 <- as.numeric(dimnames(arrays[["mastercp"]])[[2]])
  #  if (is.null(psis)) psis <- as.numeric(dimnames(arrays[["mastercp"]])[[3]])
  #  if (is.null(alph)) alph <- (1 - as.numeric(dimnames(arrays[["mastercp"]])[[5]])/100)
  #    names(longlab) <- dimnames(cis)[[3]]

  if (!is.null(psis)) par3 <- psis
  if (!is.null(phis)) par3 <- rep_len(phis, length.out = length(p1))
#  p.grid <- limits[1] + ppoints(n.grid, a = 0.5) * (limits[2] - limits[1])
#  nn <- qbinom(0.99999, n, limits[2]) # Upper limit for n when considering Poisson rates?
#  p1 <- p2 <- p.grid

#  mastercp <- array(NA, dim = c(length(p1), length(p2), length(par3), length(mymethods), length(alph), 6, 1, 1))
  mastercp <- array(NA, dim = c(length(p1), length(mymethods), length(alph), 6, 1, 1))
  dimnames(mastercp) <- list(paste(p2),
#                             paste(p2),
#                             paste(par3),
                             mymethods,
                             paste(100*(1-alph)),
                             c("cp", "lncp", "rncp", "mncp", "len", "locindex"),
                             paste(n),
                             contrast)


#  arrays <- list(xs = xs, cis = cis, mastercp = mastercp, summaries = summaries)

#  n.grid <- length(p2)
  set.seed(2012) # ensure we use the same jitters for each run
#  px <- pxg <- expand.grid(p1, p2, par3)
  px <- pxg <- cbind(p1, p2, par3)


#  refine <- 1 / n.grid
#  if (jitt) {
#    jitter <- (limits[2] - limits[1]) * array(runif(2 * dim(pxg)[1], -refine / 2, refine / 2), dim = dim(pxg)) # OK if n.grid is even
#    px[, 1:2] <- pxg[, 1:2] + jitter[, 1:2]
#  } else {
    px <- pxg
#  }

  # dimnames(px)[2]=c("p1","p2")
  if (contrast %in% c("RR")) {
    theta <- (px[, 1] / px[, 2])
  } else if (contrast == "OR") {
    #    theta <- log(px[, 1] * (1 - px[, 2]) / (px[, 2] * (1 - px[, 1])))
    theta <- (px[, 1] * (1 - px[, 2]) / (px[, 2] * (1 - px[, 1])))
  } else if (contrast == "RD") {
    theta <- px[, 1] - px[, 2]
  } # theta estimate for every p1,p2 pair

  # can we improve efficiency here? -
  # not without creating an array which is too large for available RAM.


  # alpha <- 0.05
  for (alpha in alph) {
    print(paste0("alpha=", alpha))
    #    if (contrast != "OR") {
    ci <- cis[,,, paste(100 * (1 - alpha)),,]
    #    } else ci <- log(cis[,,, paste(100 * (1 - alpha)),])

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
#i <- 1
    pb <- txtProgressBar(min = 0, max = dim(px)[1], style = 3) #text based bar
    for (i in 1:dim(px)[1]) {
      setTxtProgressBar(pb, i)
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
      # - could loop by px and then by alpha?

      # exclude (0,0) CI from CP calcs by setting its probability to zero and
      # rescaling all other probabilities - not advisable
      # prob<-prob*(1-prob[1])
      # prob[1]<-0

      #      cbind(xs, prob[, i])
      #      prob[prob == 0] <- NA

      if (any(!is.na(prob))) {

        cpl[i, ] <- t(ci[, 1, ] <= theta[i] & ci[, 2, ] >= theta[i] & ci[, 2, ] > ci[, 1, ]) %*% prob # 2-sided coverage probability. NB degenerate intervals excluded
        #        cpl[i, ] <-  t(prob) %*% (ci[, 1, ] <= theta[i] & ci[, 2, ] >= theta[i] & ci[, 2, ] > ci[, 1, ])
        lncpl[i, ] <- t(ci[, 1, ] > theta[i] | ci[, 2, ] == ci[, 1, ]) %*% prob # L-sided non-coverage (R-side is a mirror image)
        rncpl[i, ] <- t(ci[, 2, ] < theta[i] | ci[, 2, ] == ci[, 1, ]) %*% prob # R-sided non-coverage

        # not really interested in length, I think location is more important
        # - and its complicated to explain on log scale
        # but included as reviewers might request it
        lens <- ci[, 2, ] - ci[, 1, ]
        ## you get pretty whacky results with length on linear scale
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

    # Free up some memory
    rm(ci)

    # convert CPs etc to a square grid of p1,p2 points
#    mydim <- c(length(p1), length(p2), length(par3), nmeth)
#    cp <- array(cpl, dim = mydim)
#    lncp <- array(lncpl, dim = mydim)
#    rncp <- array(rncpl, dim = mydim)
#    mncp <- array(mncpl, dim = mydim)
#    len <- array(lenl, dim = mydim)
    # Calculate Newcombe's location index
#    locindex <- ifelse(cp == 1, 0, mncp / (1 - cp))
#    dimnames(cp) <- dimnames(locindex) <- dimnames(lncp) <-
#      list(paste(p1), paste(p2), paste(par3), mymethods)
##      dimnames(rncp) <- dimnames(len) <- dimnames(mncp) <-

    dimnames(cpl) <- dimnames(locindexl) <- dimnames(lncpl) <-
      dimnames(rncpl) <- dimnames(lenl) <- dimnames(mncpl) <-
      list(p2, mymethods)


#    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "cp", , ] <- cp
#    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "lncp", , ] <- lncp
#    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "rncp", , ] <- rncp
    #    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "mncp", , ] <- mncp
#    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "locindex", , ] <- locindex
#    arrays[["mastercp"]][, , , , paste(100 * (1 - alpha)), "len", , ] <- len

    mastercp <- array(rbind(cpl, lncpl, rncpl, mncpl, locindexl, lenl), dim = c(length(p2), 6, nmeth))
    dimnames(mastercp) <- list(p2, c("cp", "lncp", "rncp", "mncp", "locindex", "len") , mymethods)

    # Free up some memory
    rm(
  #    cp, lncp, rncp, mncp, locindex, len,
       cpl, lncpl, rncpl, mncpl, lenl, locindexl)

  }

  aperm(mastercp, c(1,3,2))
}


#save(cparrays_RD40, file = paste(outpath, "cparrays.RD.40.200.Rdata", sep = ""))

if (FALSE) {

# to retrieve a previously run array:
# outpath <- "/Users/ssu/Documents/Main/Courses_papers/skewscore/paired/"
# load(file=paste(outpath, "cparrays.RR.", 60, ".",100,".Rdata",sep=""))
# arr <- myarrays

#### RR
# n=20
# 4min + 17min
# 4min + 47min with 200 psps
# 3min with skewcorr +
# 1.5min for cis per alpha incl skew
# +1.5min per combo for n=100 and no smoothing
# +1.7min per combo for n=100 with smoothing
# +6min per combo with n=200 and no smoothing
# +7min per combo with n=200 with smoothing
# n=40
# 25min, +2.25h for 2x2
# 11min per alpha (33)
# +39min per combo for n=200 with smoothing (8h for 3*4)
system.time(arrays_RD40 <- cifun(n=40, contrast="RD", n.grid=200,
                                 #alph = c(0.01, 0.05, 0.1), psis = c(1, 2, 10, 100)))[[3]]/60
        #                         alph = c(0.01, 0.05, 0.1), psis = c(1, 2, 10, 100)))[[3]]/60
                                 alph = c(0.05), psis = c(10)))[[3]]/60
Sys.time()
system.time(cparrays_RD40 <- cpfun(arrays = arrays_RD40, contrast="RD", smooth=TRUE))[[3]]/60

system.time(arrays_RR40 <- cifun(n=40, contrast="RR", n.grid=200,
                                 #alph = c(0.01, 0.05, 0.1), psis = c(1, 2, 10, 100)))[[3]]/60
                                 alph = c(0.01, 0.05, 0.1), psis = c(1, 2, 10, 100)))[[3]]/60
Sys.time()
system.time(cparrays_RR40 <- cpfun(arrays = arrays_RR40, contrast="RR", smooth=TRUE))[[3]]/60

teamlist <- list(RRpairteam, RR2pairteam, RRccpairteam)
teamlabels <- c("RRpair", "RR2pair", "RRccpair")

teamlist <- list(RDpairteam, RD2pairteam, RDccpairteam)
teamlabels <- c("RDpair", "RD2pair", "RDccpair")

for (j in c(10)) {
  for (i in c(0.05)) {
    for (k in 1) {
      plotpanel(plotdata=cparrays_RD40, alpha=i, psi=j, nums="40",
                limits=c(0,1),sel=teamlist[[k]],oneside=F,plotlab=teamlabels[k],
                res.factor=6,fmt="png",linesx=F,colour=T,sided="R",
                smoothed=FALSE)
    }
  }
}


}






if (FALSE) {

# 25min + ?? for n=40
# 25min + 6h with 200 psps
# 30min with sc's + 37min with 100 psps x 2x2
system.time(arrays_RR40 <- cifun(n=40, contrast="RR", n.grid=100,
#                                 alph = c(0.01, 0.05, 0.1), psis = c(2, 10, 100)))[[3]]/60
                                 alph = c(0.05, 0.1), psis = c(2, 10)))[[3]]/60
Sys.time()
system.time(cparrays_RR40 <- cpfun(arrays = arrays_RR40, contrast="RR"))[[3]]/60
arr <- cparrays_RR40
arr[["summaries"]]["2", , "95",,]

# (was) 4 hours for N=60
#60min cis + 34min for 1 cp calc with 100 psps
#146min cis inc skewcorr + 6.5h 100grid
system.time(arrays_RR60 <- cifun(n=60, contrast="RR", n.grid=100,
                                 alph = c(0.05), psis = c(2)))[[3]]/60
#                                  alph = c(0.01, 0.05, 0.1), psis = c(1, 2, 10, 100)))[[3]]/60

Sys.time()
system.time(cparrays_RR60 <- cpfun(arrays = arrays_RR60, contrast="RR"))[[3]]/60
#system.time(cparrays_RR60 <- cpfun(arrays = cparrays_RR60, contrast="RR"))[[3]]/60
arrays <- cparrays_RR60

#### RD
# 4min + 35min without smoothing
system.time(arrays_RD20 <- cifun(n=20, contrast="RD", n.grid=100,
#                                 alph = c(0.01, 0.05, 0.1), psis = c(2, 10, 100)))[[3]]/60
                              alph = c(0.05), psis = c(2)))[[3]]/60
Sys.time()
system.time(cparrays_RD20 <- cpfun(arrays = arrays_RD20, contrast="RD"))[[3]]/60
#arr <- cparrays_RD20

# 29min + 12h for n=40
# (rollaply step responsible for long times?)
# 14min + 8min for 1 cp calc with 100grid
# 14min + 31min with 200grid
system.time(arrays_RD40 <- cifun(n=40, contrast="RD", n.grid=200,
                                 alph = 0.05, psis = c(10)))[[3]]/60
Sys.time()
system.time(cparrays_RD40 <- cpfun(arrays = arrays_RD40, contrast="RD",
                                   psis = 10, alph = 0.05))[[3]]/60
arr <- cparrays_RD40
}

if (FALSE) {
# quick test runs with n=10
system.time(arrays_RR10 <- cifun(n=10, contrast="RR", n.grid=100,
                                 alph = c(0.05, 0.1),  psis = c(2)))[[3]]/60
Sys.time()
system.time(cparrays_RR10 <- cpfun(arrays = arrays_RR10, contrast = "RR"))[[3]]/60
arr <- cparrays_RR10

system.time(arrays_RD10 <- cifun(n=10, contrast="RD", n.grid=100,
                                 alph = c(0.01, 0.05, 0.1), psis = c(2, 10, 100)))[[3]]/60
system.time(cparrays_RD10 <- cpfun(arrays = arrays_RD10, contrast="RD"))[[3]]/60
arr <- cparrays_RD10

system.time(arrays_OR10 <- cifun(n=10, contrast="OR", psis = c(2, 10, 100)))[[3]]/60
Sys.time()
system.time(cparrays_OR10 <- cpfun(arrays = arrays_OR10, contrast = "OR"))[[3]]/60
arr <- cparrays_OR10

}


if (FALSE) {
#### OR

system.time(arrays_OR20 <- cifun(n=20, contrast="OR", psis = c(2, 10, 100)))[[3]]/60
Sys.time()
system.time(cparrays_OR20 <- cpfun(arrays = arrays_OR20, contrast = "OR"))[[3]]/60
arr <- cparrays_OR20
arrays <- cparrays_OR20


system.time(arrays_OR40 <- cifun(n=40, contrast="OR", psis = c(2, 10, 100)))[[3]]/60
Sys.time()
system.time(cparrays_OR40 <- cpfun(arrays = arrays_OR40, contrast = "OR"))[[3]]/60
arr <- cparrays_OR40
}














if (FALSE) {

###############################
## Old stuff to dispose of below

par(mfrow=c(2,4), pty='s', mar=(c(2,2,2,0)+0.1))

#psi <- 2
psi <- 10
alp <- 0.05
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Tang",lside=F,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Tang-bc",lside=F,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-nw",lside=F,lines=T,CIlen=F, res.factor=1.5)
#CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-nj",lside=F,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-ns",lside=F,lines=T,CIlen=F, res.factor=1.5)
#CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="TDAS",lside=F,lines=T,CIlen=F, res.factor=1.5)

CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Tang", avg = T, lside=F,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Tang-bc", avg = T,lside=F,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-nw", avg = T,lside=F,lines=T,CIlen=F, res.factor=1.5)
#CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-nj", avg = T,lside=F,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-ns", avg = T,lside=F,lines=T,CIlen=F, res.factor=1.5)


CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Tang", lside=T,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Tang-bc", lside=T,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-nw", lside=T,lines=T,CIlen=F, res.factor=1.5)
#CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-nj", lside=T,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-ns", lside=T,lines=T,CIlen=F, res.factor=1.5)
#CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="TDAS", lside=T,lines=T,CIlen=F, res.factor=1.5)


CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Tang", lside=T, avg = T,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Tang-bc", lside=T, avg = T,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-nw", lside=T, avg = T,lines=T,CIlen=F, res.factor=1.5)
#CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-nj", lside=T,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-ns", lside=T, avg = T,lines=T,CIlen=F, res.factor=1.5)
#CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="TDAS", lside=T,lines=T,CIlen=F, res.factor=1.5)

abline(coef = c(0, 0.5))

CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Tang-ccdr",lside=F,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Tang-cc125",lside=F,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-ccnw",lside=F,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-ccns",lside=F,lines=T,CIlen=F, res.factor=1.5)

CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Tang-ccdr", lside=T,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Tang-cc125", lside=T,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-ccnw", lside=T,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="MOVER-ccns", lside=T,lines=T,CIlen=F, res.factor=1.5)



CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="SCAS",lside=F,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Jeffreys",lside=F,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="mid-p",lside=F,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Wilson",lside=F,lines=T,CIlen=F, res.factor=1.5)

CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="SCAS", lside=T,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Jeffreys", lside=T,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="mid-p", lside=T,lines=T,CIlen=F, res.factor=1.5)
CPcontour(plotdata=arr, alpha=alp, psi=psi, nums=1, methlab="Wilson", lside=T,lines=T,CIlen=F, res.factor=1.5)



CPcontour(plotdata=arr,alpha=0.05,phi=0.5,nums=1,methlab="Tango",lside=F,lines=T,CIlen=F)
CPcontour(plotdata=arr,alpha=0.05,phi=0.2,nums=1,methlab="Tang",lside=F,lines=T,CIlen=F)
CPcontour(plotdata=arr,alpha=0.05,phi=0.0,nums=1,methlab="Tang",lside=F,lines=T,CIlen=F)

system.time(arr <- runarrays(n=30,n.grid=100,phis=c(0,0.2,0.5),contrast="RD"))[3]/60
#arr
CPcontour(plotdata=arr,alpha=0.05,phi=0.0,nums=1,methlab="Tango",lside=T,lines=T,CIlen=F)
CPcontour(plotdata=arr,alpha=0.05,phi=0.5,nums=1,methlab="Tango",lside=T,lines=T,CIlen=F)
CPcontour(plotdata=arr,alpha=0.05,phi=0.2,nums=1,methlab="Tango",lside=T,lines=T,CIlen=F)

CPcontour(plotdata=arr,alpha=0.05,phi=0.,nums=1,methlab="TDAS",lside=F,lines=T,CIlen=F)
CPcontour(plotdata=arr,alpha=0.05,phi=0.2,nums=1,methlab="Tang",lside=F,lines=T,CIlen=F)


}






if (FALSE) {

another <- function(m,
                    n,
                    contrast = FALSE,
                    exact = FALSE,
                    mbc = FALSE,
                    n.grid = 200,
                    alph = c(0.01,0.05,0.1,0.2),
                    prerun = FALSE,
                    plots = TRUE,
                    lines = F,
                    res.factor = 6,
                    fmt = "png",
                    n.grid.old = 200,
                    colour = T,
                    limits = c(0,1),
                    plotlim = c(0,1),
                    sided = "R"
                    ) {

	refine <- 1/n.grid
	p.grid <- limits[1]+ppoints(n.grid)*(limits[2]-limits[1])

if (contrast %in% c("ID","IDRR")){
#	mm<-qpois(0.99999,m*limits[2])
	nn<-qpois(0.99999,n*limits[2])
} else {
#	mm<-qbinom(0.99999,m,limits[2])
	nn<-qbinom(0.99999,n,limits[2])
}

	p1 <- p2 <- p.grid #+jitter #add random jitter here if required
	if (contrast=="OR") {
	  methods<-dimnames(diffBinconf.all.OR(5,5,10,10))[[1]]
	} else if (contrast=="RR") {
	  methods<-dimnames(diffBinconf.all.RR(5,5,10,10))[[1]]
	} else if (contrast=="ID") {
	  methods<-c("MN","SC","MOVER-J","Wa","SCcc","SCcomp","MOVER-E","MOVER-comp")
	} else if (contrast=="IDRR") {
	  methods<-c("MN","SC","SCcc","SCcomp","Jef","MOVER-J","MOVER-E","MOVER-comp","Wa","BC-U","BC-J")
	} else {
	  methods<-dimnames(diffBinconf.all(5,5,10,10))[[1]]
	}

	mastercp <- array(NA,dim=c(length(p1),length(p2),length(methods),4,6,1))
	dimnames(mastercp) <- list(paste(p1),paste(p2),methods,c(80,90,95,99),c("cp","lncp","avecp","avelncp","len"))#,"pow"))
	summaries <- array(NA,dim=c(length(methods),4,22,1))
	dimnames(summaries) <- list(methods,c(80,90,95,99),c("meanlen","meanCP","minCP","pctCons","pctBad","pctnear","pctAveCons","pctAvenear","pctAvenearlo","pctAvenearhi","pctAveBad","MNCP","DNCP","maxLNCP","pctBad.1side","pctCons.1side","pctnear.1side","pctAvenear.1side","typeI","maxtypeI","avetypeI","maincolour"))
	cis <- array(NA,dim=c(length(methods),3,(mm+1)*(nn+1),4,1))
	dimnames(cis)[[4]] <- c(80,90,95,99)
	dimnames(cis)[[2]] <- c("LCL","UCL","disjoint")
	#myarrays<-list(mastercp,summaries,cis)
	#rm(mastercp)
	#rm(summaries)
	#rm(cis)

	if(prerun==TRUE) {
	  load(file=paste(outpath,"masterarraysx.",ifelse(sided=="L","L.",""),ifelse(contrast=="OR","OR.",ifelse(contrast=="RR","RR.",ifelse(contrast=="ID","ID.",ifelse(contrast=="IDRR","IDRR.","")))),m,".",n,".",n.grid,".Rdata",sep=""))
		#dimnames(myarrays[[1]])[[3]]<-dimnames(myarrays[[2]])[[1]]<-dimnames(myarrays[[3]])[[1]]<-methods
#
#		dimnames(myarrays[[3]])[[3]]<-c("meanlen","meanCP","minCP","pctCons","pctBad","pctnear","pctAveCons","pctAvenear","pctAvenearlo","pctAvenearhi","pctAveBad","MNCP","DNCP","maxLNCP","pctBad.1side","pctCons.1side","pctnear.1side","pctAvenear.1side","typeI","maxtypeI","avetypeI","maincolour")
	#	myarrays<-bigfun(m=m,n=n,p1=p1,p2=p1,alph=alph,arrays=myarrays,exact=exact,window=0.08,prerun=prerun,sided=sided)
	#	save(myarrays,file=paste(outpath,"masterarrays.",m,".",n,".Rdata",sep=""))
	} else if(prerun=="cis"|prerun=="some") {  load(file=paste(outpath,"masterarraysx.",ifelse(sided=="L","L.",""),ifelse(contrast=="OR","OR.",ifelse(contrast=="RR","RR.",ifelse(contrast=="ID","ID.",ifelse(contrast=="IDRR","IDRR.","")))),m,".",n,".",n.grid.old,".Rdata",sep=""))
#		dimnames(myarrays[[1]])[[3]]<-dimnames(myarrays[[2]])[[1]]<-dimnames(myarrays[[3]])[[1]]<-methods
		cis<-myarrays[[3]]
#		cis[1:30,,,,]<-myarrays[[3]]

		if(prerun %in% c("cis","some")) myarrays<-list(mastercp,summaries,cis)	#overwrite existing CP grid to allow a change in n.grid without recalculating all the CIs. If prerun="some", nothing is wiped (unless extra methods are being added)
		rm(mastercp); rm(summaries); rm(cis)
system.time(myarrays <- bigfun(m=m,n=n,mm=mm,nn=nn,p1=p1,p2=p1,alph=alph,arrays=myarrays,exact=exact,window=0.08,contrast=contrast,prerun=prerun,limits=limits,sided=sided)
)
		save(myarrays,file=paste(outpath,"masterarraysx.",ifelse(sided=="L","L.",""),ifelse(contrast=="OR","OR.",ifelse(contrast=="RR","RR.",ifelse(contrast=="ID","ID.",ifelse(contrast=="IDRR","IDRR.","")))),m,".",n,".",n.grid,".Rdata",sep=""))
	} else {
		myarrays<-list(mastercp,summaries,cis)
		rm(mastercp); rm(summaries); rm(cis)
		(myarrays<-bigfun(m=m,n=n,mm=mm,nn=nn,p1=p1,p2=p1,alph=alph,arrays=myarrays,exact=exact,window=0.08,contrast=contrast,limits=limits,sided=sided))
		save(myarrays,file=paste(outpath,"masterarraysx.",ifelse(sided=="L","L.",""),ifelse(contrast=="OR","OR.",ifelse(contrast=="RR","RR.",ifelse(contrast=="ID","ID.",ifelse(contrast=="IDRR","IDRR.","")))),m,".",n,".",n.grid,".Rdata",sep=""))
	}

#  n <- 40
#	nums=paste(m,",",n,sep="")
	nums <- paste(n)
	RRpairteam <- c("Tang", "MOVER-w", "MOVER-nw","MOVER-ns") 	#Paired RR

	ateam <- c("MN","Mee","Wbc","SC") 	#Score 4
	bteam <- c("Wald","AC","BL","N") 	#Simple 4
	cteam <- c("JP","AMB","MOVER-J","NJa") #"Jeff") #Jeffreys-based
	dteam <- c("GNcc","Ncc","HA","Waldcc") #continuity corrected 4
	eteam <- c("AME","CZ","SS","NE") 		#exact 4
	fteam <- c("LR","GN","MNBL","MNBLcc")  #some others
	mpteam <- c("AMmp","CZmp","SSmp","NEmp")  #some others

#	mbcteam<-c("MN","MN2","MNBL","Wbc2") #experimental modified bias correction (shelved)
#	mnblteam<-c("MN","MNBL","MNBL60","MNBL55") #different weightings for MNBL
	bigteam<-c("SC","MN","BL","N") #for large sample section of paper
	pageteam<-c("N","MN","AC","MNBL") #for 1-page document


	#for skewscore paper:
	skewteam<-c("SC","MN","MOVER-J","Wald")
	IDteam<-c("SC","MN","MOVER-J","Wa")
	RRteam<-c("SC","MN","MOVER-J","Katz") #RR scale intervals (with contrast="RR")
	RRIDteam<-c("SC","MN","MOVER-J","Wa")
	#ORteam<-c("SC","MN","MOVER-J","Gart") #logit scale intervals (with contrast="OR)
	ORteam<-c("SC","MN","MOVER-J","Woolf") #logit scale intervals (with contrast="OR)

	#for supplement section comparing Bayesian methods:
	#(RD versions to be added)
	RDteam2<-c("MOVER-J","AMB","AMA","NM")
	RDteam3<-c("SCcc","SCcomp","MOVER-E","MOVER-comp")
	RRteam2<-c("SC","MOVER-J","AMB","NM")
	RRteam3<-c("SCcc","SCcomp","MOVER-E","MOVER-comp")
	RRteam4<-c("SC","SC","Katz","Gart")
	RRIDteam2<-c("SC","MOVER-J","BC-J","BC-U")
	RRIDteam3<-c("SCcc","SCcomp","MOVER-E","MOVER-comp")
	IDteam3<-c("SCcc","SCcomp","MOVER-E","MOVER-comp")
	ORteam2<-c("SC","MOVER-J","AMB","NM")
	ORteam3<-c("SCcc","SCcomp","MOVER-E","MOVER-comp")
	ORteam4<-c("SCnobcf","Agresti","Woolf","Gart")

#	plotpanel(plotdata=myarrays,alpha=0.05,nums=nums,sel=pageteam,oneside=F,plotlab="1page",res.factor=res.factor,fmt="png")
#	plotpanel(plotdata=myarrays,alpha=0.05,nums=nums,sel=pageteam,oneside=T,plotlab="1page",res.factor=res.factor,fmt="png")
if (plots == TRUE) {
#alpha=0.05
myarrays <- arrays_RR40

	for(alpha in alph) {
		if (contrast=="OR") {
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=ORteam,oneside=F,plotlab="OR",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=ORteam,oneside=T,plotlab="OR",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=ORteam2,oneside=F,plotlab="OR2",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=ORteam2,oneside=T,plotlab="OR2",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=ORteam4,oneside=F,plotlab="OR4",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=ORteam4,oneside=T,plotlab="OR4",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=ORteam3,oneside=F,plotlab="OR3",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=ORteam3,oneside=T,plotlab="OR3",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
		} else if(contrast=="RR"){

		  plotpanel(plotdata=cparrays_RR40,alpha=0.05, psi=10,nums=1,limits= plotlim,sel=RRpairteam,oneside=F,plotlab="RRpair",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
		  plotpanel(plotdata=cparrays_RR40,alpha=0.05, psi=10,nums=1,limits= plotlim,sel=RRpairteam,oneside=T,plotlab="RRpair",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)

		  plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=RRteam,oneside=F,plotlab="RR",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
		  plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=RRteam,oneside=T,plotlab="RR",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=RRteam2,oneside=F,plotlab="RR2",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=RRteam2,oneside=T,plotlab="RR2",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=RRteam4,oneside=F,plotlab="RR4",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=RRteam4,oneside=T,plotlab="RR4",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=RRteam3,oneside=F,plotlab="RR3",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=RRteam3,oneside=T,plotlab="RR3",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)

#plotpanel(CIlen=T,plotdata=myarrays,alpha=alpha,nums=nums,limits= plotlim,sel=RRteam,oneside=F,plotlab="RRlen",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour)

		} else if(contrast=="ID"){
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=IDteam,oneside=F,plotlab="ID",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=IDteam,oneside=T,plotlab="ID",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=IDteam3,oneside=F,plotlab="ID3",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=IDteam3,oneside=T,plotlab="ID3",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
		} else if(contrast=="IDRR"){
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=RRIDteam,oneside=F,plotlab="IDRR",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=RRIDteam,oneside=T,plotlab="IDRR",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=RRIDteam2,oneside=F,plotlab="IDRR2",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=RRIDteam2,oneside=T,plotlab="IDRR2",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=RRIDteam3,oneside=F,plotlab="IDRR3",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=RRIDteam3,oneside=T,plotlab="IDRR3",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
		} else {
#			if((nums=="205,205" | nums=="315,105") & fmt=="tiff")	{
#				plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=bigteam,oneside=F,plotlab="big",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#				plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=bigteam,oneside=T,plotlab="big",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#			}
#			if(fmt=="tiff")	{
				plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=skewteam,oneside=F,plotlab="RD",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
				plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=skewteam,oneside=T,plotlab="RD",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
				plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=RDteam2,oneside=F,plotlab="RD2",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
				plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=RDteam2,oneside=T,plotlab="RD2",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
				plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=RDteam3,oneside=F,plotlab="RD3",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
				plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=RDteam3,oneside=T,plotlab="RD3",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#			}
#			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=ateam,oneside=F,plotlab="a",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel= ateam,oneside=T,plotlab="a",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#if(fmt=="tiff") {
#plotpanel(CIlen=T,plotdata=myarrays,alpha=alpha,nums=nums,sel=ateam,oneside=F,plotlab="alen",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#plotpanel(CIlen=T,plotdata=myarrays,alpha=alpha,nums=nums,sel=bteam,oneside=F,plotlab="blen",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#}
#			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=bteam,oneside=F,plotlab="b",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=bteam,oneside=T,plotlab="b",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=cteam,oneside=F,plotlab="c",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=cteam,oneside=T,plotlab="c",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=dteam,oneside=F,plotlab="d",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=dteam,oneside=T,plotlab="d",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)

#plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=fteam,oneside=F,plotlab="f",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=fteam,oneside=T,plotlab="f",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			if(exact==TRUE) {
				plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=eteam,oneside=F,plotlab="e",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
				plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,limits=plotlim,sel=eteam,oneside=T,plotlab="e",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#
#				#plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=mpteam,oneside=F,plotlab="mp",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#				#plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=mpteam,oneside=T,plotlab="mp",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			}
##			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=eteam,oneside=F,plotlab="e",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
##			plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=eteam,oneside=T,plotlab="e",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#			if(mbc==TRUE) {
#				plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=mbcteam,oneside=F,plotlab="mbc",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#				plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=mbcteam,oneside=T,plotlab="mbc",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
#			}
			#plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=mnblteam,oneside=F,plotlab="m",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
			#plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=mnblteam,oneside=T,plotlab="m",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
		}
	}
}

#output indication of disjoint intervals for exact methods
if(exact) myarrays[[3]][c("AME"),3,,,][apply(myarrays[[3]][c("AME"),3,,,],1,sum)>0,]

}

}

if(FALSE) {



alpha=c(0.01,0.05,0.1)
#alpha <- 0.05
m=30;n=30
m=45;n=15
#m=15;n=45

#need to run these again for IDRR
m=150;n=50
m=100;n=100
m=50;n=150

m=13;n=27
m=20;n=10
m=7;n=3
m=5;n=5
m=15;n=15
m=205;n=205
m=105;n=315
m=315;n=105

#plots for skewscore paper
outpath=paste(root,"Main/Courses_papers/skewscore/plots/",sep="")
system.time(another(m,n,logit=F,exact=F,lines=F,alph=alpha,prerun=T,fmt="png"))[[3]]/60 #
system.time(another(m,n,logit="ID",exact=F,limits=c(0,1.05),plotlim=c(0,1),n.grid=210,lines=F,alph=alpha,prerun=T,fmt="png"))[[3]]/60 #
system.time(another(m,n,logit="RR",exact=F,lines=F,alph=alpha,prerun=F,fmt="tiff"))[[3]]/60 #
system.time(another(m,n,logit="IDRR",exact=F,limits=c(0,1.05),plotlim=c(0,1),n.grid=210,lines=F,alph=alpha,prerun=T,fmt="png"))[[3]]/60 #
system.time(another(m,n,logit="OR",exact=F,lines=F,alph=alpha,prerun=T,fmt="tiff"))[[3]]/60 #


#combine plots for publication
#install.packages("png")
library(png)
fig1a <- readPNG(paste0(outpath,"_png/summaryRD95_150,50os.png"), TRUE)
fig1b <- readPNG(paste0(outpath,"_png/summaryID95_150,50os.png"), TRUE)
res.factor <- 6
tiff(file=paste0(outpath,"_tiff/Figure1.tiff"),width=600*res.factor,height=2*330*res.factor,type="quartz")
#png(file=paste0(outpath,"_png/Figure1.png"),width=600*res.factor,height=2*330*res.factor,type="quartz")
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(0:1,0:1, type='n',axes=F, xaxs="i", yaxs="i")
#plot(0:1,0:1, axes=F, xaxs="i", yaxs="i")
rasterImage(fig1a, 0, 0.5, 1, 1)
rasterImage(fig1b, 0, 0, 1, 0.5)
text(0.02,0.98,"(a)",cex=1.5*res.factor)
text(0.02,0.48,"(b)",cex=1.5*res.factor)
dev.off()

fig2a <- readPNG(paste0(outpath,"_png/summaryRR95_150,50os.png"), TRUE)
fig2b <- readPNG(paste0(outpath,"_png/summaryIDRR95_150,50os.png"), TRUE)
res.factor <- 6
tiff(file=paste0(outpath,"_tiff/Figure2.tiff"),width=600*res.factor,height=2*330*res.factor)
#png(file=paste0(outpath,"_png/Figure2.png"),width=600*res.factor,height=2*330*res.factor)
par(mar=c(0,0,0,0))
plot(0:1,0:1, type='n',axes=F, xaxs="i", yaxs="i")
rasterImage(fig2a, 0, 0.5, 1, 1)
rasterImage(fig2b, 0, 0, 1, 0.5)
text(0.02,0.98,"(a)",cex=1.5*res.factor)
text(0.02,0.48,"(b)",cex=1.5*res.factor)
dev.off()



#fig4 <- readPNG("/Users/petelaud/Documents/Main/Courses_papers/skewscore/Laud2015/Figure4.eps", TRUE)


#zoomed plots for skewscore paper  limits=c(0,0.05),plotlim=c(0,0.05),
#m=150;n=50
m=600;n=200
#m=500;n=2500
outpath=paste(root,"Main/Courses_papers/skewscore/plots/zoom/",sep="")
system.time(another(m,n,logit=F,exact=F,lines=F,alph=alpha,prerun=F,fmt="tiff",limits=c(0,0.05),plotlim=c(0,0.05)))[[3]]/60 #
system.time(another(m,n,sided="L",logit=F,exact=F,lines=F,alph=alpha,prerun=T,fmt="tiff",limits=c(0,0.05),plotlim=c(0,0.05)))[[3]]/60 #
system.time(another(m,n,logit="OR",exact=F,lines=F,alph=alpha,prerun=F,fmt="tiff",limits=c(0,0.05),plotlim=c(0,0.05)))[[3]]/60 #
system.time(another(m,n,logit="OR",sided="L",exact=F,lines=F,alph=alpha,prerun=T,fmt="tiff",limits=c(0,0.05),plotlim=c(0,0.05)))[[3]]/60 #
system.time(another(m,n,logit="RR",exact=F,lines=F,alph=alpha,prerun=F,fmt="tiff",limits=c(0,0.05),plotlim=c(0,0.05)))[[3]]/60 #
system.time(another(m,n,logit="ID",exact=F,n.grid=210,lines=F,alph=alpha,prerun=F,fmt="tiff",limits=c(0,0.0525),plotlim=c(0,0.05)))[[3]]/60 #
system.time(another(m,n,logit="IDRR",exact=F,n.grid=210,lines=F,alph=alpha,prerun=F,fmt="tiff",limits=c(0,0.0525),plotlim=c(0,0.05)))[[3]]/60 #

load(file=paste(outpath,"masterarraysx.IDRR.",m,".",n,".",210,".Rdata",sep=""))
load(file=paste(outpath,"masterarraysx.RR.",m,".",n,".",200,".Rdata",sep=""))
myarrays[[3]][,,,"95",]
myarrays[[2]][,"95",,]
dimnames(myarrays[[2]])
min(myarrays[[1]][,,5,"95","lncp",])

ratesCI(6,3,11000,6000,dist="poi",contrast="RR")




### Tables for publication
outpath=paste(root,"Main/Courses_papers/skewscore/plots/",sep="")
outorder<-c("SC","MN","MOVER-J","Wald")
#numslist<-c("15,15","30,30","45,45","20,10","45,15","60,30","205,205","173,346","315,105") #)
numslist<-c("30,30","45,15","100,100","150,50","50,150")
m=n=30
n.grid=200
load(file=paste(outpath,"masterarraysx.",m,".",n,".",n.grid,".Rdata",sep=""))
#methods<-outorder
summaries<-array(NA,dim=c(length(outorder),5,2,length(numslist),5)) #[c(2,1,3,4)])
#methods<-dimnames(diffBinconf.all(5,5,10,10))[[1]]
dimnames(summaries)[[1]]<-outorder
dimnames(summaries)[[2]]<- c("bRD","pRD","bRR","pRR","bOR")  #dimnames(myarrays[[2]])[2:3]
dimnames(summaries)[[4]]<-numslist
dimnames(summaries)[[5]] <- c("pctnear","pctAvenear","pctnear.1side","pctAvenear.1side","pctAveCons")
lev=c("90","95")
dimnames(summaries)[[3]]<-lev
for(nums in numslist) {
	mn<-strsplit(nums,",")[[1]]
	n.grid=200
	load(file=paste(outpath,"masterarraysx.",mn[1],".",mn[2],".",n.grid,".Rdata",sep=""))
	summaries[outorder,"bRD",,nums,]<-myarrays[[2]][outorder,lev,c("pctnear","pctAvenear","pctnear.1side","pctAvenear.1side","pctAveCons"),]
	load(file=paste(outpath,"masterarraysx.RR.",mn[1],".",mn[2],".",n.grid,".Rdata",sep=""))
	outorder2<-c("SC","MN","MOVER-J","Katz")
	summaries[outorder,"bRR",,nums,]<-myarrays[[2]][outorder2,lev,c("pctnear","pctAvenear","pctnear.1side","pctAvenear.1side","pctAveCons"),]
	load(file=paste(outpath,"masterarraysx.OR.",mn[1],".",mn[2],".",n.grid,".Rdata",sep=""))
#	outorder2<-c("SC","MN","MOVER-J","Gart")
	outorder2<-c("SC","MN","MOVER-J","Woolf")
	summaries[outorder,"bOR",,nums,]<-myarrays[[2]][outorder2,lev,c("pctnear","pctAvenear","pctnear.1side","pctAvenear.1side","pctAveCons"),]
	n.grid=210
	load(file=paste(outpath,"masterarraysx.ID.",mn[1],".",mn[2],".",n.grid,".Rdata",sep=""))
	outorder2<-c("SC","MN","MOVER-J","Wa")
	summaries[outorder,"pRD",,nums,]<-myarrays[[2]][outorder2,lev,c("pctnear","pctAvenear","pctnear.1side","pctAvenear.1side","pctAveCons"),]
	load(file=paste(outpath,"masterarraysx.IDRR.",mn[1],".",mn[2],".",n.grid,".Rdata",sep=""))
	summaries[outorder,"pRR",,nums,]<-myarrays[[2]][outorder2,lev,c("pctnear","pctAvenear","pctnear.1side","pctAvenear.1side","pctAveCons"),]
}

summaries[,4,,,"pctAvenear.1side"]
summaries[,,3:6,"pctnear.1side"]

library(xtable)

xtable(((summaries[,1,"95",,c("pctAvenear.1side")])))
xtable(((summaries[,2,"95",,c("pctAvenear.1side")])))
xtable(((summaries[,3,"95",,c("pctAvenear.1side")])))
xtable(((summaries[,4,"95",,c("pctAvenear.1side")])))
xtable(((summaries[,5,"95",,c("pctAvenear.1side")])))


ftable(noquote(summaries[,,5:6,c("pctAvenear.1side")]),col.vars=c(3,2))



ftable(noquote(summaries[,,3:6,c("pctnear.1side")]),col.vars=c(3,2))
#Table 1
ftable(noquote(summaries[,,c("95","90"),c("pctAvenear","pctnear","pctAvenear.1side","pctnear.1side")]),col.vars=c(3,2))






dev.off()
prerun<-"some"
prerun<-F
exact<-F
prerun<-T
exact<-T

#diffBinconf.all(2,3,7,3,exact=F)
###Plots for website
									#timings with exact=F:
t1<-system.time(x1<-another(5,5,exact=exact,prerun=prerun)); t1/60; range(c(x1)[abs(c(x1))>0])  #2min / 24min(exact) / 6.8(prerun) / 18.8(prerun=cis) #39min (searchnum=30) no disjoint intervals
#17min

t2<-system.time(x2<-another(7,3,exact=exact,prerun=prerun)); t2/60; range(c(x2)[abs(c(x2))>0])  #2min / 26min(exact) / (prerun)  disjoint range: 0.036,0.063 / 0.073,0.107 as proportion of CI width
#15min /11(exact=F)
t3<-system.time(x3<-another(15,15,exact=exact,prerun=prerun,res.factor=6)); t3/60; range(c(x3)[(c(x3))>0 & abs(c(x3))<1]) #5min  / 6.2(prerun=T) / 122min(exact=T) /103  /213 searchnum=30 disjoint range: 0.015, 0.122
#80min / 5(prerun) /19(exact=F)
#t3<-system.time(x3<-another(15,14,exact=exact,prerun=prerun,res.factor=6)); t3/60; range(c(x3)[(c(x3))>0 & abs(c(x3))<1]) #5min  / 6.2(prerun=T) / 122min(exact=T) /103  /213 searchnum=30 disjoint range: 0.015, 0.122
t4<-system.time(x4<-another(20,10,exact=exact,prerun=prerun)); t4/60; range(c(x4)[abs(c(x4))>0 & abs(c(x4))<1])  # 114min(exact) / 6min(prerun)  /81min(searchnum=10) 140min(searchnum=30) disjoint range: 0.014, 0.178
#90min /21(exact=F)

#prerun<-T
#exact<-F
t5<-system.time(x5<-another(30,30,exact=exact,prerun=prerun)); t5/60; range(c(x5)[abs(c(x5))>0 & abs(c(x5))<1])  #69min 6.2(prerun) 19.5h(exact) 32h searchnum=30  max disjoint: 0.016
#27h /44min(exact=F)
t6<-system.time(x6<-another(45,15,exact=exact,prerun=prerun)); t6/60; range(c(x6)[abs(c(x6))>0 & abs(c(x6))<1])  #77min 6.2(prerun) 11h+ (exact) 13h searchnum=30  disjoint range: 0.015,0.225 (obtained with wide=0.2)
#look at the distribution of x6 to see if disjoint intervals tend to be shorter or longer
#12h /42min(exact=F)
#7.8h


system.time(another(45,45,exact=F,prerun=prerun))/60  #166min
#93min (exact=F) 3.3(prerun)
system.time(another(60,30,exact=F,prerun=prerun))/60  #160min
#91min (exact=F) 3.3(prerun)

#system.time(another(7,3,alph=0.1,exact=T,prerun=F,n.grid=200))[[3]]/60  #2min / 4.5min(300) / 28min(exact) / 6.8(prerun)
#system.time(another(30,30,alph=c(0.2),exact=T,prerun=F,n.grid=200))[[3]]/60  #2min / 4.5min(300) / 28min(exact) / 6.8(prerun)

#re-run with searchnum=20:
#system.time(another(45,15,alph=c(0.2),exact=T,prerun=F,n.grid=200))[[3]]/60  #2min / 4.5min(300) / 28min(exact) / 6.8(prerun)
#system.time(another(30,30,alph=c(0.05),exact=T,prerun=F,n.grid=200))[[3]]/60  #2min / 4.5min(300) / 28min(exact) / 6.8(prerun)

####Colour & B&W Plots for publication
system.time(another(15,15,exact=T,prerun=T,lines="2",fmt="tiff",alph=0.05,colour=T))[[3]]/60 # / 1.5(prerun=T) / 23min(exact=T) #222 searchnum=30
system.time(another(20,10,exact=T,prerun=T,lines="2",fmt="tiff",alph=0.05,colour=T))[[3]]/60  #21min(exact)

#timer<-system.time(another(173,346,prerun=F,alph=c(0.05),n.grid=200)); timer[[3]]/60  #469
#450min
#timer<-system.time(another(205,205,exact=F,prerun=prerun,n.grid=200)); timer[[3]]/60  #112min n.grid=100 /
#31h(prerun=cis)
#timer<-system.time(another(315,105,exact=F,prerun=prerun,n.grid=200)); timer[[3]]/60  #112min n.grid=100 / 470min n.grid=200 /2.5min xx(prerun=100)  4x6.5h or maybe less?(prerun=cis)
#18h
#system.time(another(315,105,prerun=T,n.grid=100))[[3]]/60  #4x262min n.grid=200 / 15.5h n.grid=100



##########################
#Colour Plots for publication
system.time(another(205,205,prerun=T,alph=0.05,lines="2",fmt="tiff",colour=T))[[3]]/60  #112min n.grid=100 / 470min n.grid=200
#timer<-system.time(another(173,346,prerun=F,alph=c(0.05),plots=F,n.grid=200)); timer[[3]]/60  #469
#450min
system.time(another(315,105,prerun=T,alph=0.05,lines="2",fmt="tiff",colour=T))[[3]]/60  #262min n.grid=200
#system.time(another(105,315,prerun=F,lines="both",mbc=F,alph=c(0.05),n.grid=200,fmt="tiff"))[[3]]/60  #827min? with power off n.grid=200

system.time(another(346,173,logit="RR",limits=c(0,1),exact=F,lines=F,alph=0.05,prerun=F,fmt="tiff"))[[3]]/60 #
system.time(another(346,173,logit="ID",limits=c(0,0.5),exact=F,lines=F,alph=0.05,prerun=F,fmt="tiff"))[[3]]/60 #
system.time(another(346,173,logit="IDRR",limits=c(0,0.5),exact=F,lines=F,alph=0.05,prerun=F,fmt="tiff"))[[3]]/60 #

system.time(another(105,1005,logit="ID",limits=c(0,0.1),exact=F,lines=F,alph=0.05,prerun=F,fmt="tiff"))[[3]]/60 #
system.time(another(5000,5000,logit="ID",limits=c(0,0.01),exact=F,lines=F,alph=0.05,prerun=F,fmt="tiff"))[[3]]/60 #
system.time(another(500,500,logit="ID",limits=c(0,0.1),exact=F,lines=F,alph=0.05,prerun=F,fmt="tiff"))[[3]]/60 #
system.time(another(50,50,logit="ID",limits=c(0,1),exact=F,lines=F,alph=0.05,prerun=F,fmt="tiff"))[[3]]/60 #
system.time(another(300,600,logit="ID",limits=c(0,0.1),exact=F,lines=F,alph=0.05,prerun=F,fmt="tiff"))[[3]]/60 #

system.time(another(50,50,logit="IDRR",limits=c(0,1),exact=F,lines=F,alph=0.05,prerun=F,fmt="tiff"))[[3]]/60 #






system.time(another(1000,1000,logit="ID",limits=c(0,0.015),plotlim=c(0,0.01),exact=F,lines=F,alph=0.05,prerun=F,fmt="tiff"))[[3]]/60 #
dev.off()

system.time(another(15,15,logit="RR",limits=c(0,1),exact=F,lines=F,alph=0.05,fmt="tiff"))[[3]]/60 #

system.time(another(15,15,logit=F,limits=c(0,0.1),exact=F,lines=F,alph=0.05,prerun="cis",fmt="tiff"))[[3]]/60 #

jkl<-(diffBinconf.all.RR(1,4,100,100))
(jkl[,2,])/(jkl[,1,])

diffBinconf.RRID(1,0,50,50,skew=F)
1/diffBinconf.RRID(0,1,50,50,skew=F)
diffBinconf.RRID(1,0,50,50,skew=T)
1/diffBinconf.RRID(0,1,50,50,skew=T)
diffBinconf.RRID(c(0,1,0),c(0,0,1),500,500)

exp(diffBinconf.all.RR(5,50,500,500))


max(rbinom(100000,100,0.1)) #need to restrict the CIs being calculated when limits are (0,0.1)

#system.time(another(15,15,logit=T,exact=F,,lines="2",alph=0.05,prerun=F,fmt="tiff"))[[3]]/60 #
#system.time(another(30,15,logit=T,exact=F,lines="2",n.grid=200,prerun=F,mbc=F))[[3]]/60 #5min / 10.2(300) / 3.6(prerun=T) / 102min(exact=T)
#system.time(another(30,30,logit=T,exact=F,lines="2",n.grid=200,prerun=F,mbc=F))[[3]]/60 #5min / 10.2(300) / 3.6(prerun=T) / 102min(exact=T)
#system.time(another(101,101,logit=T,exact=F,lines="2",n.grid=100,alph=0.05,prerun=F,mbc=F))[[3]]/60 #


#system.time(another(40,40,alph=0.05,prerun=F,lines="2"))[[3]]/60  #107min
#system.time(another(50,50,prerun=F,mbc=T))[[3]]/60  #40min / 85
#system.time(another(60,20,prerun=F,n.grid=100))[[3]]/60  #9min / 19
#system.time(another(90,45,prerun=F))[[3]]/60   #139min
#system.time(another(90,90,prerun=F,mbc=T))[[3]]/60   #123 / 271	 #127min
#system.time(another(120,60,prerun=T,mbc=T))[[3]]/60   #113min / 243
#system.time(another(80,120))[[3]]/60		 #150min / 320

#system.time(another(162,108,prerun=F,alph=0.05,n.grid=200,fmt="tiff"))[[3]]/60  #262min n.grid=200
#dimnames(myarrays[[1]])[[3]]<-dimnames(myarrays[[2]])[[1]]<-methods

#locate minCP
m=100
n=100
n.grid=210
load(file=paste(outpath,"masterarraysx.IDRR.",m,".",n,".",n.grid,".Rdata",sep=""))
load(file=paste(outpath,"masterarrays.",m,".",n,".",n.grid,"alph5.Rdata",sep=""))
myarrays[[3]]["SC",1:2,,"95",]
#[myarrays[[3]]["NE","disjoint",,"95",]>0]
ratesCI(20,0,10,10,contrast="RR",dist="poi")
check1<-myarrays[[1]][,,"MN","95","cp",]
check1<-myarrays[[1]][,,"Wa","95","cp",];
#dimnames(
exp(myarrays[[3]]["MN",,,"95",])
sum(is.na(myarrays[[1]][,,"AM",,"lncp",]))
myarrays[[2]][,"95",,]
sum(is.na(myarrays[[3]][,"LCL",,"90",]))

myarrays[[3]]["NE","disjoint",,"95",][myarrays[[3]]["NE","disjoint",,"95",]>0]

selrow<-apply(check1,1,function(x) sum(x<0.852))
selcol<-apply(check1,2,function(x) sum(x<0.852))
check1[names(selrow)[selrow>0],names(selcol)[selcol>0]]
check1[names(selcol)[selcol>0],names(selrow)[selrow>0]]
check1[1:10,200-(0:9)]  #minimum is not right in the corner, but very close to it.
check1[200-(0:9),(1:10)]


#Match Gart Nam numbers
m=315; n=105
m=n=205
n.grid=200
load(file=paste(outpath,"masterarrays.",m,".",n,".",n.grid,".Rdata",sep=""))

myarrays[[1]][,,c("MN","SC"),"95","lncp",]
myarrays[[3]]["MN",1:2,1:3,"95",]
myarrays[[3]]["SC",1:2,1:3,"95",]
diffBinconf.RRID(1,0,50,50,skew=F)
myarrays[[1]][c("0.0125","0.2025"),c("0.2025","0.6025"),c("MN","GN"),"95","lncp",]
myarrays[[1]][c("0.2025","0.6025"),c("0.0125","0.2025"),c("MN","GN"),"95","lncp",]

#Data supporting statement in paper about MN not controlling Type I error
round(unlist(myarrays[[1]]["0.8025","0.9025",c("MN","GNbc","MNBL","LR","Wald","AC"),"95","lncp",]),5)

p1<-as.numeric(dimnames(myarrays[[1]][,,c("MN"),"95","avelncp",])[[1]])
p2diag<-paste(p1+0.1)[p1+0.1<=1]
mths<-dimnames(myarrays[[1]])[[3]]
typeIerr<-sapply(mths,function(i) diag(myarrays[[1]][paste(p1),p2diag,i,"95","avelncp",]))

plot(p1[p1+0.1<=1],typeIerr[,"BL"],type="l",ylim=c(0.02,0.03))
abline(h=0.025*(1+c(-0.1,0,0.1)))
for (i in 1:(length(mths)-11)) {
#i<-0
#i<-i+1
lines(p1[p1+0.1<=1],typeIerr[,mths[-c(3,22,24,27,29,17,18,19,20,21,30)][i]])
#mths[i]
}
plot(p1[p1+0.1<=1],typeIerr[,"MN"],type="l",ylim=c(0.02,0.03))
abline(h=0.025*(1+c(-0.1,0,0.1)))
for (i in 1:5) {
#i<-0
#i<-i+1
lines(p1[p1+0.1<=1],typeIerr[,mths[c(3,22,24,27,29)][i]])
#mths[c(3,22,24,27,29)][i]

}

#for(m)


check30<-cbind(myarrays[[3]]["NE","LCL",,"95",]-myarrays[[3]]["CZ","LCL",,"95",],expand.grid(0:30,0:30),t(myarrays[[3]]["CZ",c("LCL","disjoint"),,"95",]))

check30<-cbind(myarrays[[3]]["CZ","UCL",,"95",]-myarrays[[3]]["CZ","LCL",,"95",],expand.grid(0:30,0:30))
plot(check30[,2]-check30[,3],check30[,1])

pos30<-check30[check30[,1]< -0.015,]
contour(0:30,0:30,matrix(pos30[,1],c(31,31)))
diffBinconf.exact(18,4,30,30,method="NE")
diffBinconf.exact(18,4,30,30,method="CZ")

#looking into odd zones of low LNCP for CZ
myarrays[[1]]["0.5025","0.3025","CZ","95","avelncp",]
lls<-myarrays[[3]]["CZ","LCL",,"95",]
lls[abs(lls-0.2)<0.01]



### output table
numslist<-c("15,15","30,30","45,45","20,10","45,15","60,30")#,"205,205","173,345","315,105")
outorder=c("MN","BL","MNBL","Wbc","AC","N","Mee","Wald","MNBL55","MNBL","MNBL66")
#, "20,40","50,50","90,90","120,60","80,120")
#numslist<-c("205,205","173,345","315,105")
#outorder=c("MN","Mee","BL","MNBL","Wbc","W","AC","N","PL","Wald") #,"MN2","Wbc2","MN2BL")

m=n=50
n.grid=200
load(file=paste(outpath,"masterarrays.IDRR.",m,".",n,".",n.grid,".Rdata",sep=""))
load(file=paste(outpath,"masterarrays.",30,".",30,".",200,".Rdata",sep=""))
load(file=paste(outpath,"masterarrays.",45,".",15,".",n.grid,".Rdata",sep=""))
summaries<-array(NA,dim=c(dim(myarrays[[2]])[1:3],length(numslist))[c(1,4,2,3)])
dimnames(summaries)[c(1,3,4)]<-dimnames(myarrays[[2]])[1:3]
dimnames(summaries)[[2]]<-numslist
for(nums in numslist) {
	mn<-strsplit(nums,",")[[1]]
	load(file=paste(outpath,"masterarrays.",mn[1],".",mn[2],".",n.grid,".Rdata",sep=""))
	summaries[outorder,nums,,]<-myarrays[[2]][outorder,,,]
}
summaries.small<-summaries


dimnames(myarrays[[2]])[[1]]
dimnames(myarrays[[1]])[[3]]
dimnames(myarrays[[3]])[[1]]

load(file=paste(outpath,"masterarrays.",30,".",30,".",200,".exact.Rdata",sep=""))
myarrays[[2]]["NE","95",,]
load(file=paste(outpath,"masterarrays.",30,".",30,".",200,".bak.Rdata",sep=""))
myarrays[[2]]["NE","95",,]

myarrays[[3]]["NE","LCL",,"95",]

### Tables for publication
outorder<-c("Wald","AC","BL","N","MN","Mee","Wbc","GNbc","LR","MNBL")
#outorder=c("MN","BL","MNBL","GN","GNbc","Wbc","AC","N","Mee","LR","Wald")
#outorder=c("MN","BL","MNBL","Wbc","AC","N","Mee","LR","Wald","JP","AMB","NJ")
#(,"MNCC","MNBLCC","NCC","HA","WaldCC","AME","CZ","SS","NE")
numslist<-c("15,15","30,30","45,45","20,10","45,15","60,30","205,205","173,346","315,105") #)
#numslist<-c("205,205","173,345","315,105")
m=n=15
n.grid=200
load(file=paste(outpath,"masterarrays.",m,".",n,".",n.grid,".Rdata",sep=""))
#methods<-outorder
summaries<-array(NA,dim=c(length(outorder),dim(myarrays[[2]])[2:3],length(numslist))[c(1,4,2,3)])
#methods<-dimnames(diffBinconf.all(5,5,10,10))[[1]]
dimnames(summaries)[[1]]<-outorder
dimnames(summaries)[c(3,4)]<-dimnames(myarrays[[2]])[2:3]
dimnames(summaries)[[2]]<-numslist
for(nums in numslist) {
	mn<-strsplit(nums,",")[[1]]
	load(file=paste(outpath,"masterarrays.",mn[1],".",mn[2],".",n.grid,".Rdata",sep=""))
#	dimnames(myarrays[[2]])[[1]]<-methods
	summaries[outorder,nums,,]<-myarrays[[2]][outorder,,,]
}
#summaries.large<-summaries


#Add brackets to avenear summaries indicating anticonservative methods
anticons<-as.numeric(summaries[,,,"pctAveCons"])<50
anticons2<-as.numeric(summaries[,,,"pctCons"])<50
anticons[is.na(anticons)]<-FALSE
anticons2[is.na(anticons2)]<-FALSE
summaries[,,,"pctAvenear"]<-paste(ifelse(anticons,"("," "),summaries[,,,"pctAvenear"],ifelse(anticons,")"," "),sep="")
summaries[,,,"pctnear"]<-paste(ifelse(anticons2,"("," "),summaries[,,,"pctnear"],ifelse(anticons2,")"," "),sep="")

save(summaries,file=paste(outpath,"publishedsummaries.Rdata",sep=""))



#Table 1
ftable(noquote(summaries[,1:6,c("95","90"),c("pctAvenear")]),col.vars=c(3,2))
#Table 2
ftable(noquote(summaries[,1:6,c("95","90"),c("pctAvenear.1side")]),col.vars=c(3,2))
#Table 1
ftable(noquote(summaries[,1:6,c("95","90"),c("pctAvenear","pctnear","pctAvenear.1side","pctnear.1side")]),col.vars=c(3,2))

anticons==anticons2


#confs<-c("99","95","90","80")
#confs<-c("95","90")
#sink(file=paste(outpath,"summary2012f.txt",sep=""))
#,outpath,summaries,outorder
#dimnames(summaries)
outtable<-function(conf){
	sink(file=paste(outpath,"summary",".txt",sep=""),type="output")
	cat("Table 1: pctAvenear\n")
	print(ftable(noquote(summaries[outorder,1:6,conf,"pctAvenear"]),col.vars=c(3,2)))
	cat("\n\n")
	cat("Table 2: pctAvenear.1side\n")
	print(ftable(noquote(summaries[outorder,1:6,conf,"pctAvenear.1side"]),col.vars=c(3,2)))
	cat("\n\n")
	cat("Table 3: large samples\n")
	print(ftable(noquote(summaries[outorder,7:9,"95",c("pctAvenear","pctAvenear.1side")]),col.vars=c(3,2)))
	cat("\n\n")
	cat("Table S1: pctnear\n")
	print(ftable(noquote(summaries[outorder,1:6,conf,"pctnear"]),col.vars=c(3,2)))
	cat("\n\n")
	cat("Table S2: pctnear.1side\n")
	print(ftable(noquote(summaries[outorder,1:6,conf,"pctnear.1side"]),col.vars=c(3,2)))
	cat("\n\n")
	cat("Table S3: large samples, alpha=0.05\n")
	print(ftable(noquote(summaries[outorder,7:9,"95",c("pctnear","pctnear.1side")]),col.vars=c(3,2)))
	cat("\n\n")
	cat("Table S4: Additional large sample summaries, alpha=0.05\n")
	print(ftable(noquote(summaries[outorder,7:9,"95",c("pctBad.1side","minCP")]),col.vars=c(3,2)))
	cat("\n\n")
#	print("pctAveCons")
#	print(ftable(summaries[outorder,,conf,"pctAveCons"],col.vars=c(2)))
	sink()
}
outtable(conf=c("95","90"))


outtable("80")
outtable("90")
outtable("99")



#confs<-c("99","80")
#sink(file=paste(outpath,"summary2012b.txt",sep=""))
#"pctAvecons"
#summaries[outorder,,c("95","90"),"maincolour"]




"pctAvenear"
ftable(summaries[outorder,,confs,"pctAvenear"],col.vars=c(3,2))
"pctAvecons"
ftable(summaries[outorder,,confs,"pctAveCons"],col.vars=c(3,2))
"pctnear"
ftable(summaries[outorder,,confs,"pctnear"],col.vars=c(3,2))
"pctBad.1side"
ftable(summaries[outorder,,confs,"pctBad.1side"],col.vars=c(3,2))
"pctAvenear.1side"
ftable(summaries[outorder,,confs,"pctAvenear.1side"],col.vars=c(3,2))
"pctnear.1side"
ftable(summaries[outorder,,confs,"pctnear.1side"],col.vars=c(3,2))
"minCP"
ftable(summaries[outorder,,confs,"minCP"],col.vars=c(3,2))
sink()


conf="90"
conf="95"
alph=1-as.numeric(conf)/100
nums<-"15,15"
nums<-"20,10"
minx<-0; pos=4
nums<-"205,205"
nums<-"315,105"
minx<-50; pos=2


outorder=c("MN","BL","MNBL","Wbc","AC","N","LR","Wald","GNbc")

#Fig 6b plotting MACP % proximate
sumryplot<-function(conf,nums,minx,pos){
	alph=1-as.numeric(conf)/100
	x<-as.numeric(gsub("[^0-9.]","",summaries[outorder,nums,conf,"pctAvenear"]))
	y<-as.numeric(gsub("[^0-9.]","",summaries[outorder,nums,conf,"pctAvenear.1side"]))
	z<-as.numeric(summaries[outorder,nums,conf,"pctAveCons"])<50
	plot(x,y,xlim=c(minx,100),ylim=c(0,100),main=paste("N=(",nums,"), alpha=",alph,sep=""),xaxs="i",yaxs="i",xlab="",ylab="",pty="s")
	text(x,y,labels=paste(ifelse(z,"(",""),outorder,ifelse(z,")",""),sep=""),pos=pos,offset=0.25)
	mtext(side=1,text="MACP % proximate",line=2,cex=res.factor*0.8)
	mtext(side=2,text="1-sided MACP % proximate",line=2,cex=res.factor*0.8)
	abline(h=seq(0,100,10),v=seq(0,100,10),lty=3,col="gray",xpd=F)
}
res.factor=2
tiff(file=paste(outpath,"fig6b.GN.tif",sep=""),width=400*res.factor,height=600*res.factor,compression="none")
par(mfrow=c(3,2))
par(mar=c(4,3,2,1)+0.1,xpd=T,cex=1.5,pty='s',tcl=-0.4)
sumryplot("95","15,15",0,c(4,2,4,4,4,2,4,4,4))
sumryplot("95","20,10",0,c(4,2,4,4,2,4,4,4,4))
sumryplot("90","15,15",0,c(4,4,4,4,4,4,4,4,3))
sumryplot("90","20,10",0,4)
sumryplot("95","205,205",50,c(2,2,2,2,4,2,4,2,3))
sumryplot("95","315,105",50,c(4,4,2,2,2,2,4,2,3))
dev.off()

#Alternative version plotting %proximate
sumryplot<-function(conf,nums,minx,pos){
	alph=1-as.numeric(conf)/100
	x<-as.numeric(gsub("[^0-9.]","",summaries[outorder,nums,conf,"pctnear"]))
	y<-as.numeric(gsub("[^0-9.]","",summaries[outorder,nums,conf,"pctnear.1side"]))
	z<-as.numeric(summaries[outorder,nums,conf,"pctCons"])<50
	plot(x,y,xlim=c(minx,100),ylim=c(0,100),main=paste("N=(",nums,"), alpha=",alph,sep=""),xaxs="i",yaxs="i",pty="s",xlab="",ylab="")
#	text(x,y,labels=outorder,pos=pos)
	text(x,y,labels=paste(ifelse(z,"(",""),outorder,ifelse(z,")",""),sep=""),pos=pos,offset=0.25)
	mtext(side=1,text="% proximate",line=2,cex=res.factor*0.8)
	mtext(side=2,text="1-sided % proximate",line=2,cex=res.factor*0.8)
	abline(h=seq(0,100,10),v=seq(0,100,10),lty=3,col="gray",xpd=F)
}
res.factor=2
tiff(file=paste(outpath,"fig6bx.GN.tif",sep=""),width=400*res.factor,height=600*res.factor,compression="none")
par(mfrow=c(3,2))
par(mar=c(4,3,2,1)+0.1,xpd=T,cex=1.5,pty='s',tcl=-0.4)
sumryplot("95","15,15",0,c(2,4,4,4,2,4,4,4,3))
sumryplot("95","20,10",0,c(4,2,4,4,2,4,4,4,3))
sumryplot("90","15,15",0,c(2,4,4,4,4,4,2,4,3))
sumryplot("90","20,10",0,c(4,4,4,4,4,4,4,4,3))
sumryplot("95","205,205",50,c(2,2,4,3,4,2,2,2,3))
sumryplot("95","315,105",50,c(2,4,4,2,2,2,2,2,3))
dev.off()





dimnames(myarrays[[2]])
myarrays[[2]][c("MN","N","MN/BLJ","AM"),"95","minCP",] #minCP
myarrays[[2]][c("MN","N","MN/BLJ","AM"),"95","pctBad",] #minCP

nums=paste(m,",",n,sep="")
oneside=F
alpha=0.05
plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("N","MN","BLJ","MN/BLJ"),oneside=oneside,plotlab="a")

myarrays[[3]][,,,paste(100*(1-alpha)),]


dimnames(myarrays[[1]])[[3]]<-dimnames(myarrays[[2]])[[1]]<-methods

myarrays[[1]][1,1,1,1,1,1]

#identify which objects are taking up most memory
#objects.length <- function(list)
#{
#   res <- sapply(list, function(x) round(object.size(get(x))/1000000))
#   names(res) <- list
#   res
#}
#objects.length(ls())


#rm(myarrays)

##n.grid=200 #199 #odd number of grid points produces more symmetrical plots
##refine=1/n.grid
##limits=c(0,1)
##xlim=ylim=limits #Choose range of p1,p2 to be studied - for 'zooming in' on the plot, e.g. xlim=ylim=c(0.7,1)
##avedelta=0.05 #used later for averaged plots
##xlim.ext=c(max(0+refine/2,(xlim[1]-avedelta)),min(1-refine/2,(xlim[2]+avedelta))) #shrink limits by refine/2 to get midpoints of grid
##ylim.ext=c(max(0+refine/2,(ylim[1]-avedelta)),min(1-refine/2,(ylim[2]+avedelta))) #extended limits to allow calculation of mean CP for p+/-avedelta (only applies if limits are not (0,1))
##p1.grid=seq(xlim.ext[1],xlim.ext[2],refine)
##jitter1=runif(n.grid,-refine/2,refine/2)
##jitter2=runif(n.grid,-refine/2,refine/2)
###jitter=c(jitter1,rev(jitter1))
##jitter1=jitter2=0 #add zero jitter to get precise values of delta
##p1=p1.grid + jitter1 #adding random jitter avoids bias in summary stats
##p2.grid=seq(ylim.ext[1],ylim.ext[2],refine)
##p2=p2.grid + jitter2
####p1=sort(runif(150,min=xlim.ext[1],max=xlim.ext[2])) #alternative random uniform sampling produces ugly plots
####p2=sort(runif(150,min=ylim.ext[1],max=ylim.ext[2]))


##methods=c("MN","Wbc","AC","N","Mee","Wall","MN2","Wbc2","Wald","PL","NCC","WaldCC","AM","CZ","EE","JP","BL")#
#methods=c("OR1","OR2","OR3","MN")
#rm(mastercp)
#mastercp.Unif=array(NA,dim=c(n.rand,16,6,4,3))
#dimnames(mastercp.Unif)=list(1:n.rand,methods,1:6,c(80,90,95,99),c("cp","lncp","len"))
#dimnames(mastercp.Unif)[[3]]<-c("5,5","15,15","25,35","20,50","50,50","90,90")

#dimnames(summaries)=list(methods,c("5,5","15,15","25,35","20,50","50,50","90,90","80,120"),c(80,90,95,99),c("meanlen","meanCP","pctCons","pctBad","pctnear","pctAveCons","pctAvenear","MNCP","DNCP","pctnear.1side","pctAvenear.1side"))
#dimnames(summaries)=list(methods,c("5,5","15,15","30,20","20,40","90,90","120,60","80,120"),c(80,90,95,99),c("meanlen","meanCP","pctCons","pctBad","pctnear","pctAveCons","pctAvenear","MNCP","DNCP","pctnear.1side","pctAvenear.1side","typeI","avetypeI"))
#summaries[,,,1:12]<-myarrays[[2]]
#myarrays<-list(myarrays[[1]],summaries)

#mastercp[,,,1:6,,]<-myarrays[[1]][,,,1:6,,]
#summaries[,1:6,,]<-myarrays[[3]]


system.time(myarrays<-bigfun(m=15,n=15,alpha=alpha,p1=p1,p2=p1,arrays=myarrays,exact=F))[3]/60 #2min
save(myarrays,file=paste(outpath,"masterarrays.15.15.Rdata",sep=""))



#for cross-checking against 2-d plots in papers, e.g. Brown&Li, Fagerland:
m=40;n=15
p.grid<-ppoints(200)
p.grid<-seq(0,1,0.01)
#p2<-c(0.2,0.5) #horizontal slice
p2<-p.grid #diagonal slice
p1<-p.grid #add random jitter here if required
methods<-dimnames(diffBinconf.all(5,5,10,10))[[1]]
mastercp=array(NA,dim=c(length(p1),length(p2),length(methods),4,6,1))
dimnames(mastercp)<-list(paste(p1),paste(p2),methods,c(80,90,95,99),c("cp","lncp","avecp","avelncp","len"))#,"pow"))
summaries<-array(NA,dim=c(length(methods),4,21,1))
#dimnames(summaries)=list(methods,c(80,90,95,99),c("meanlen","meanCP","pctCons","pctBad","pctnear","pctAveCons","pctAvenear","MNCP","DNCP","pctnear.1side","pctAvenear.1side","typeI","avetypeI"))
dimnames(summaries)=list(methods,c(80,90,95,99),c("meanlen","meanCP","minCP","pctCons","pctBad","pctnear","pctAveCons","pctAvenear","pctAvenearlo","pctAvenearhi","pctAveBad","MNCP","DNCP","maxLNCP","pctBad.1side","pctnear.1side","pctAvenear.1side","typeI","maxtypeI","avetypeI","maincolour"))
cis<-array(NA,dim=c(length(methods),2,(m+1)*(n+1),4,1))
dimnames(cis)[[4]]<-c(80,90,95,99)
myarrays<-list(mastercp,summaries,cis)
rm(mastercp)
#rm(mastercp.Unif)
rm(summaries)
system.time(slice<-bigfun(m=m,n=n,p1=p1,p2=p2,alph=0.05,arrays=myarrays,square=FALSE,jitt=FALSE))[3]/60

#Brown&Li plot showing errors in their fig 4a.
png(file=paste(outpath,"BrownLiFig4a.png",sep=""))
#plot(p1-0.5,slice[[1]][,"0.5","BLR","95","cp",1],ylim=c(0.88,0.98),type="l",ylab="CP",main=expression(paste("Coverage probability of Recentered CI for p1=0.5, m=n=20",alpha)))
#x=2
#mtext(side=3,line=1,text=substitute(paste("test",alpha,x),list(x=x)))
#mtext(side=3,line=2,text="test")

abline(h=0.95)
#lines(p1-0.5,slice[[1]][,"0.5","N","95","cp",1],lty=2)
#lines(0.5-p2,slice[[1]]["0.5",,"BLj","95","cp",1],lty=3)
dev.off()

dimnames(slice[[1]])

#Fagerland plot & why their conclusion on MN is wrong.  They've used Mee for starters.
del<-0
fmethod<-"BL"
attrib<-"cp"
#attrib<-"len"
#attrib<-"pow"
p1diag<-paste(p1[p2+del>=0 & p2+del<=1])
p2diag<-paste(p2[p1-del>=0 & p1-del<=1])
diagplot<-diag(slice[[1]][p1diag,p2diag,fmethod,"95",attrib,1])
length(diagplot)
plot(p2diag,diagplot,type="l",ylim=c(0.92,1),ylab="power",xlab="p2",main=paste("Coverage probability for NJ,\n at increments of delta from 0 to 0.3. (n1=n2=40)")) #,
#plot(p2diag,diagplot,type="l",ylim=c(0,0.05),ylab="power",xlab="p2",main=paste("Coverage probability for NJ,\n at increments of delta from 0 to 0.3. (n1=n2=40)")) #,
#plot(p1diag,diagplot,type="l",ylim=c(0.95,0.98),ylab="power") #,
abline(h=0.95)
abline(h=0.945,lty=2)
dels<-del+seq(0,-0.3,-0.01)
i=28
for(i in 1:length(dels)){
#dels[29]
	deli<-dels[i]
	p1diag<-paste(p1[round(p1-deli,3)>=0 & round(p1-deli,3)<=1])
	p2diag<-paste(p2[round(p2+deli,3)>=0 & round(p2+deli,3)<=1])
	diagplot<-diag(slice[[1]][p1diag,p2diag,fmethod,"95",attrib,1])
	lines(p2diag,diagplot,type="l",lty=1)
}

#main=expression(paste("Power of Brown-Li CI for p1=p2, m=n=20")))
lines(p1,diag(slice[[1]][,,"MN","95",attrib,1]),lty=2)
lines(p1,diag(slice[[1]][,,"Mee","95",attrib,1]),lty=2)
lines(p1,diag(slice[[1]][,,"MNBL","95",attrib,1]),lty=3)
lines(p1,diag(slice[[1]][,,"N","95",attrib,1]),lty=4)
lines(p1,diag(slice[[1]][,,"Wald","95",attrib,1]),lty=4)



#outpath="//emea.astrazeneca.net/uk/Alderley Park/Users 09/kbpj139/Documents/useful/DiffBin/plots/"
outpath="/Users/Pete/Main/Contract/ssu2010_003/DiffBin/plots/"
load(file=paste(outpath,"masterarrays.205.205.200.Rdata",sep=""))
#load(file=paste(outpath,"masterarrays.315.105.200.Rdata",sep=""))
p1=as.numeric(dimnames(myarrays[[1]])[[1]])
p2=as.numeric(dimnames(myarrays[[1]])[[2]])
methods=dimnames(myarrays[[1]])[[3]]
refine=1/length(p1)
#dimnames(myarrays[[1]])[[4]]=dimnames(myarrays[[2]])[[2]]=c("5,5","15,15","30,20","20,40","50,50","120,60","80,120")
dimnames(myarrays[[1]])[[4]]=dimnames(myarrays[[2]])[[2]]=c("16,15","31,21","41,21","51,49","91,89","89,87","121,81")
#dimnames(myarrays[[1]])[[1]]=p1


#### KEEP ####
#2-d plot of lncp at delta=-0.1
mth="LR"
plot(p1[1:180],diag(myarrays[[1]][1:180,21:200,mth,"95","avelncp",]),type="l",ylim=c(0,0.03))
abline(h=0.025)
abline(h=0.025*c(0.9,1.1),lty=2)
lines(p1[1:180],diag(myarrays[[1]][1:180,21:200,mth,"95","lncp",]),type="l",lty=2)


mydiag=diag(myarrays[[1]][1:180,21:200,"MN",nums,a,b])
maadiag=rollapply(zoo(mydiag),21,mean)
plot(p1[1:180],mydiag,type="l")
lines(p1[11:170],maadiag)

a="95"
b="lncp"
nums="91,89"
method="MN2"
plot(p1[1:180],diag(myarrays[[1]][1:180,21:200,method,nums,a,b]),type="l",ylim=c(0.015,0.045))
#lines(p1[1:180],diag(myarrays[[1]][1:180,21:200,method,"89,87",a,b]),type="l")
#lines(p1[1:180],diag(myarrays[[1]][1:180,21:200,method,"91,89",a,"avelncp"]),type="l")
abline(h=c(0.05,0.0250))
abline(h=c(0.055,0.0275),lty=2)
#lines(p1[1:180],diag(myarrays[[1]][1:180,21:200,"MN","51,49",a,"avelncp"]),lty=2)
lines(p1[1:180],diag(myarrays[[1]][1:180,21:200,"MN2",nums,a,b]),type="l",lty=2)
lines(p1[1:180],diag(myarrays[[1]][1:180,21:200,"N",nums,a,b]),type="l",lty=2)
lines(p1[1:180],diag(myarrays[[1]][1:180,21:200,"PL",nums,a,b]),type="l",lty=2)
lines(p1[1:180],diag(myarrays[[1]][1:180,21:200,"BL",nums,a,b]),type="l",lty=2)
lines(p1[1:180],diag(myarrays[[1]][1:180,21:200,"Wbc",nums,a,b]),type="l",lty=2)
lines(p1[1:180],diag(myarrays[[1]][1:180,21:200,"Wald",nums,a,b]),type="l",lty=2)


lines(p1[1:180],diag(myarrays[[1]][1:180,21:200,"PL","121,59","95",b]),type="l",lty=2)
lines(p1[1:180],diag(myarrays[[1]][1:180,21:200,"MN","121,81","95",b]),type="l",lty=2)
lines(p1[1:179],diag(myarrays[[1]][1:179,22:200,"MN","41,21","95",b]),type="l",lty=2)
lines(p1[1:181],diag(myarrays[[1]][1:181,20:200,"MN","41,21","95",b]),type="l")
lines(p1[1:178],diag(myarrays[[1]][1:178,23:200,"MN","41,21","95",b]),type="l")
lines(p1[1:182],diag(myarrays[[1]][1:182,19:200,"MN","41,21","95",b]),type="l")
lines(p1[1:180],diag(myarrays[[1]][1:180,21:200,"PL","41,21","95",b]),type="l",lty=3)
abline(h=c(0.05,0.025))
abline(h=0.0275,lty=2)
plot(p1[1:200],diag(myarrays[[1]][1:200,1:200,"MN","15,15","90",b]),type="l")
(myarrays[[1]][1:10,1:10,"Wald",2,3,"avelncp"])
#dim(myarrays[[2]])

p1[1:180]-p2[21:200]

#m=n=5; alpha=0.05
#myarrays<-bigfun(m=m,n=n,alpha=0.05,refine=refine,p1=p1,p2=p2)
#mastercp[,,,paste(m,n,sep=","),paste(100*(1-alpha)),]<-myarrays[[1]]
#mastercp.Unif[,,paste(m,n,sep=","),paste(100*(1-alpha)),]<-myarrays[[2]]

(myarrays[[1]][100,100,"EE","15,15","95","cp"])
dimnames(myarrays[[2]])

objects.length(ls())
warnings()
alpha=0.01
alpha=0.05
system.time(myarrays<-bigfun(m=5,n=5,alpha=alpha,p1=p1,p2=p2,arrays=myarrays,exact=F))/60 #68sec with n.grid=199
system.time(myarrays<-bigfun(m=15,n=15,alpha=alpha,p1=p1,p2=p1,arrays=myarrays,exact=F)) #138 sec
system.time(myarrays<-bigfun(m=30,n=20,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #6min, 8800 with exact methods
system.time(myarrays<-bigfun(m=20,n=40,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #7min
system.time(myarrays<-bigfun(m=50,n=50,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #14min
system.time(myarrays<-bigfun(m=120,n=60,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #45min
system.time(myarrays<-bigfun(m=80,n=120,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #68min
alpha=0.1
system.time(myarrays<-bigfun(m=5,n=5,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #68sec with n.grid=199
system.time(myarrays<-bigfun(m=15,n=15,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #138 sec
system.time(myarrays<-bigfun(m=30,n=20,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #344 sec, 8800 with exact methods
system.time(myarrays<-bigfun(m=20,n=40,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #383 sec
system.time(myarrays<-bigfun(m=50,n=50,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #1465 sec
system.time(myarrays<-bigfun(m=120,n=60,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #2941s (49min)
system.time(myarrays<-bigfun(m=80,n=120,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #3891s

alpha=0.05
system.time(myarrays<-bigfun(m=16,n=15,alpha=alpha,p1=p1,p2=p1,arrays=myarrays,exact=F))/60 #138 sec
system.time(myarrays<-bigfun(m=31,n=21,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #6min, 8800 with exact methods
system.time(myarrays<-bigfun(m=41,n=21,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #7min
system.time(myarrays<-bigfun(m=51,n=49,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #14min
system.time(myarrays<-bigfun(m=91,n=89,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #14min
system.time(myarrays<-bigfun(m=89,n=87,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #14min
system.time(myarrays<-bigfun(m=121,n=59,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #45min
system.time(myarrays<-bigfun(m=121,n=81,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #68min

alpha=0.1
system.time(myarrays<-bigfun(m=16,n=15,alpha=alpha,p1=p1,p2=p1,arrays=myarrays,exact=F))/60 #138 sec
system.time(myarrays<-bigfun(m=31,n=21,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #6min, 8800 with exact methods
system.time(myarrays<-bigfun(m=41,n=21,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #7min
system.time(myarrays<-bigfun(m=51,n=49,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #14min
system.time(myarrays<-bigfun(m=91,n=89,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #14min
system.time(myarrays<-bigfun(m=121,n=59,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #45min
system.time(myarrays<-bigfun(m=121,n=81,alpha=alpha,p1=p1,p2=p2,arrays=myarrays))/60 #68min


alpha=0.2
system.time(myarrays<-bigfun(m=90,n=90,alpha=alpha,p1=p1,p2=p2,arrays=myarrays)) #2941s (49min)
system.time(myarrays<-bigfun(m=80,n=120,alpha=alpha,p1=p1,p2=p2,arrays=myarrays)) #2941s (49min)

#m=n=50; alpha=0.05
#system.time(myarrays.i<-bigfun(m=m,n=n,alpha=alpha,p1=p1,p2=p2,arrays=myarrays)) #886s
#myarrays[[1]][,,,paste(m,n,sep=","),paste(100*(1-alpha)),]<-myarrays.i[[1]]
#myarrays[[2]][,paste(m,n,sep=","),paste(100*(1-alpha)),]<-myarrays.i[[2]]

#myarrays[[2]][,,"90",c("meanCP","pctnear","pctAveCons","pctAvenear","pctAvenear.1side")]
#myarrays[[2]][,"5,5","95",c("meanCP","pctnear","pctAveCons","pctAvenear","pctAvenear.1side")]

sink(file=paste(outpath,"summary.txt",sep=""))
outorder=c("Mee","MN","MN2","Wall","Wbc","Wbc2","BL","AC","N","Wald","JP","PL","NCC","WaldCC")
myarrays[[2]][outorder,,c("90","95"),c("pctAvenear","pctAveCons","pctCons")]
myarrays[[2]][outorder,,c("90","95"),"pctAvenear.1side"]
sink()

dimnames(myarrays[[1]])

save(myarrays,file=paste(outpath,"masterarrays.200grid.Rdata",sep=""))






#coarser grid for larger sample sizes:

n.grid=99 #odd number of grid points produces more symmetrical plots (???)
refine=1/n.grid
limits=c(0,1)
xlim=ylim=limits #Choose range of p1,p2 to be studied - for 'zooming in' on the plot, e.g. xlim=ylim=c(0.7,1)
avedelta=0.05 #used later for averaged plots
xlim.ext=c(max(0+refine/2,(xlim[1]-avedelta)),min(1-refine/2,(xlim[2]+avedelta))) #shrink limits by refine/2 to get midpoints of grid
ylim.ext=c(max(0+refine/2,(ylim[1]-avedelta)),min(1-refine/2,(ylim[2]+avedelta))) #extended limits to allow calculation of mean CP for p+/-avedelta (only applies if limits are not (0,1))
p1.grid=seq(xlim.ext[1],xlim.ext[2],refine)
p1=p1.grid + runif(round(1/refine),-refine/2,refine/2); #adding random jitter avoids bias in summary stats
p2.grid=seq(ylim.ext[1],ylim.ext[2],refine)
p2=p2.grid + runif(round(1/refine),-refine/2,refine/2);
####p1=sort(runif(150,min=xlim.ext[1],max=xlim.ext[2])) #alternative random uniform sampling produces ugly plots
####p2=sort(runif(150,min=ylim.ext[1],max=ylim.ext[2]))

methods=c("MN","WL","AC","N","Mee","Wall","MN*","WL*","Wald","PL","NCC","WaldCC","OR1","OR2","OR3","OR4")
rm(mastercp)
mastercp=array(NA,dim=c(length(p1),length(p2),length(methods),7,4,5))
dimnames(mastercp)=list(paste(p1),paste(p2),methods,1:7,c(80,90,95,99),c("cp","lncp","avecp","avelncp","len"))
dimnames(mastercp)[[4]]<-c("5,5","15,15","25,35","20,50","50,50","90,90","200,200")

summaries<-array(NA,dim=c(length(methods),7,4,11))
dimnames(summaries)=list(methods,c("5,5","15,15","25,35","20,50","50,50","90,90","200,200"),c(80,90,95,99),c("meanlen","meanCP","pctCons","pctBad","pctnear","pctAveCons","pctAvenear","MNCP","DNCP","pctnear.1side","pctAvenear.1side"))

myarrays.bign<-list(mastercp,summaries)
dimnames(myarrays.bign[[1]])[[4]]<-c("5,5","15,15","25,35","20,50","150,150","90,90","200,200")
dimnames(myarrays.bign[[2]])[[2]]<-c("5,5","15,15","25,35","20,50","150,150","90,90","200,200")

load(file=paste(outpath,"masterarrays.bign.Rdata",sep=""))
p1=as.numeric(dimnames(myarrays.bign[[1]])[[1]])
p2=as.numeric(dimnames(myarrays.bign[[1]])[[2]])
methods=dimnames(myarrays.bign[[1]])[[3]]


alpha=0.1
system.time(myarrays.bign<-bigfun(m=5,n=5,alpha=alpha,p1=p1,p2=p2,arrays=myarrays.bign)) #4h12min with n.grid=99
system.time(myarrays.bign<-bigfun(m=15,n=15,alpha=alpha,p1=p1,p2=p2,arrays=myarrays.bign)) #4h12min with n.grid=99
system.time(myarrays.bign<-bigfun(m=25,n=35,alpha=alpha,p1=p1,p2=p2,arrays=myarrays.bign)) #4h12min with n.grid=99
system.time(myarrays.bign<-bigfun(m=20,n=50,alpha=alpha,p1=p1,p2=p2,arrays=myarrays.bign)) #4h12min with n.grid=99
system.time(myarrays.bign<-bigfun(m=150,n=150,alpha=alpha,p1=p1,p2=p2,arrays=myarrays.bign)) #4h12min with n.grid=99
system.time(myarrays.bign<-bigfun(m=200,n=200,alpha=alpha,p1=p1,p2=p2,arrays=myarrays.bign)) #4h12min with n.grid=99

save(myarrays.bign,file=paste(outpath,"masterarrays.bign.Rdata",sep=""))

myarrays.bign[[2]][,,c("90","95"),"pctAvenear.1side"]



#this way was no quicker
#cpfun=function(px,g,meth,theta) {
#	prob=(pmf(px[1],px[2],g[,1],g[,2],m,n))
#	cp=mysum(prob*(ci[,1,meth]<=theta & ci[,2,meth]>=theta))
#	cp
#}
#cpfun(px=px[22,],g=g,meth=meth,theta=theta[22])
#sapply(1:dim(px)[1],function(i) cpfun(px=px[i,],g=g,meth=meth,theta=theta[i]),px=px,g=g,meth=meth)



### Reproducing Newcombe's method of generating PSPs, for comparable summary stats
#if (showNewc==TRUE) {
#	psi.samp=runif(10000)
#	theta.samp=runif(10000)*(1-abs(2*psi.samp-1))
#	p1.samp=psi.samp+theta.samp/2
#	p2.samp=psi.samp-theta.samp/2
#	cp.Newc=array(NA,dim=c(length(theta.samp),nmeth))
#	for (i in 1:length(theta.samp)) {
#		prob=(pmf(p1=p1.samp[i],p2=p2.samp[i],x1=g[,1],x2=g[,2],n1=m,n2=n))
#		for (meth in 1:nmeth) {
#			cp.Newc[i,meth]=mysum(prob*(ci[,1,meth]<=theta.samp[i] & ci[,2,meth]>=theta.samp[i] & ci[,2,meth]>ci[,1,meth]))
#			lncp.Newc[i,meth]=mysum(prob*(ci[,1,meth]>theta.samp[i] | ci[,2,meth]==ci[,1,meth] ))
#		}
#	}
#
#	dimnames(cp.Newc)[2]=list(methods)
#	meanCP.Newc=format(round(apply(cp.Newc,2,mean),3),nsmall=3)
#	minCP.Newc=format(round(apply(cp.Newc,2,min),3),nsmall=3)
#	pctCons.Newc=100*round(apply(cp.Newc,2,FUN=function(x) mean(x>(1-alpha))),3)
#	MNCP.Newc=format(1-round(apply(lncp.Newc[theta.samp>0,],2,mean),3),nsmall=3)
#}


#### Producing coverage summaries using random Uniform sampling for p1,p2
#p1.samp=runif(n.rand,min=xlim[1],max=xlim[2])
#p2.samp=runif(n.rand,min=ylim[1],max=ylim[2])
#theta.samp=p1.samp-p2.samp #theta estimate for every p1,p2 pair
##lor.samp=log(p1.samp*(1-p2.samp)/(p2.samp*(1-p1.samp))) #log odds ratio estimate
#cp.Unif=lncp.Unif=len.Unif=array(NA,dim=c(length(theta.samp),nmeth))
#for (i in 1:length(theta.samp)) {
#	prob=(pmf(p1=p1.samp[i],p2=p2.samp[i],x1=g[,1],x2=g[,2],n1=m,n2=n))
#	for (meth in 1:nmeth) {
#		cp.Unif[i,meth]=mysum(prob*(ci[,1,meth]<=theta.samp[i] & ci[,2,meth]>=theta.samp[i] & ci[,2,meth]>ci[,1,meth])) #2-sided coverage probability
#		lncp.Unif[i,meth]=mysum(prob*(ci[,1,meth]>theta.samp[i] | ci[,2,meth]==ci[,1,meth])) #L-sided non-coverage (R-side is a mirror image)
#		len.Unif[i,meth]=mysum(prob*(ci[,2,meth]-ci[,1,meth]))#
##		if(logit)	{
##			cp.Unif[i,meth]=mysum(prob*(ci.logit[,1,meth]<=lor.samp[i] & ci.logit[,2,meth]>=lor.samp[i] & ci.logit[,2,meth]>ci.logit[,1,meth]))
##			lncp.Unif[i,meth]=mysum(prob*(ci.logit[,1,meth]>lor.samp[i] | ci.logit[,2,meth]==ci.logit[,1,meth] ))
##			len.Unif[i,meth]=mysum(prob*(ci.logit[,2,meth]-ci.logit[,1,meth]))
##		}
#	}
#}
#mastercp.Unif.i=array(NA,dim=c(10000,16,5))
#dimnames(mastercp.Unif.i)[[3]]=c("cp","lncp","avecp","avelncp","len")
#mastercp.Unif.i[,,"cp"]<-cp.Unif
#mastercp.Unif.i[,,"lncp"]<-lncp.Unif
##mastercp.Unif[,,paste(m,n,sep=","),paste(100*(1-alpha)),"avecp"]<-avecp
##mastercp.Unif[,,paste(m,n,sep=","),paste(100*(1-alpha)),"avelncp"]<-avelncp
#mastercp.Unif.i[,,"len"]<-len.Unif
#arrays[[2]][,,paste(m,n,sep=","),paste(100*(1-alpha)),"cp"]<-cp.Unif
#arrays[[2]][,,paste(m,n,sep=","),paste(100*(1-alpha)),"lncp"]<-lncp.Unif
#arrays[[2]][,,paste(m,n,sep=","),paste(100*(1-alpha)),"len"]<-len.Unif
#meanlen=format(round(apply(len.Unif,2,mean),3),nsmall=3)
#meanCP=format(round(apply(cp.Unif,2,mean),3),nsmall=3)
#pctCons=format(100*round(apply(cp.Unif,2,FUN=function(x) mean(x>(1-alpha))),3),nsmall=1)
#pctnear=format(100*round(apply(cp.Unif,2,FUN=function(x) mean(x>(1-alpha)-0.005 & x<(1-alpha)+0.005)),3),nsmall=1)
#MNCP=format(round(apply(lncp.Unif,2,function(x) mean(x[theta.samp>0])),3),nsmall=3)
#DNCP=format(round(apply(lncp.Unif,2,function(x) mean(x[theta.samp<0])),3),nsmall=3)
##pctCons.1side=format(100*round(apply(lncp.Unif,2,FUN=function(x) mean(x<(alpha/2))),3),nsmall=1)
#pctnear.1side=format(100*round(apply(lncp.Unif,2,FUN=function(x) mean(x>(alpha/2)-0.005 & x<(alpha/2)+0.005)),3),nsmall=1)


}
