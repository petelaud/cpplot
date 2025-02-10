set.seed(2012) #ensure we use the same jitters for each run



if (FALSE) {

  # The following code reproduces the results in the manuscript
  # 'Equal-tailed confidence intervals for paired binomial proportions' by Peter J. Laud
#  pak::pak("petelaud/ratesci-dev") # Will need updating to the public repository
#  library(ratesci)

  # Set path for output files as required
#  outpath <- '/'

  root <- "/Users/ssu/Documents/"
  outpath <- paste(root, "Main/Courses_papers/skewscore/paired/", sep = "")


  ### OPTIONAL:
  ### Run the CP calculation function for N=20, N=40 and N=65
  ### WARNING: for N=65, these take several hours to run!
  RDpairteam <- RRpairteam <- c("SCAS-bc", "SCAS", "AS", "MOVER-NJ", "MOVER-W", "BP")
  alphas <- c(0.1, 0.05, 0.01)
  phis <- c(0.1, 0.25, 0.5, 0.75)
  system.time(mycis <- cifun(n=20, contrast="RD", alph = alphas))[[3]]/60
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
  system.time(mycis <- cifun(n=20, contrast="RR", alph = alphas))[[3]]/60
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
  system.time(mycis <- cifun(n=40, contrast="RD", alph = alphas))[[3]]/60
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
  system.time(mycis <- cifun(n=40, contrast="RR", alph = alphas))[[3]]/60
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
  system.time(mycis <- cifun(n=65, contrast="RD", alph = alphas, methods = RDpairteam))[[3]]/60
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
  system.time(mycis <- cifun(n=65, contrast="RR", alph = alphas, methods = RRpairteam))[[3]]/60
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60

  # Combine output arrays for summarising across different Ns
  mynums <- c(20, 40, 65)
  mymethods <- c("SCAS-bc", "SCAS", "AS", "MOVER-NJ", "MOVER-W", "BP")
  nmeth <- length(mymethods)
  load(file=paste0('data/', "cparrays.RD.", 40, ".",200,".Rdata"))
  mydims <- dim(arrays$summaries)
  mydims[5] <- length(mynums)
  mydims[6] <- 2
  mydims[2] <- nmeth
  mydimnames <- dimnames(arrays$summaries)
  mydimnames[[5]] <- paste(mynums)
  mydimnames[[6]] <- c("RD", "RR")
  mydimnames[[2]] <- mymethods

  bigarray <- array(NA, dim = mydims)
  dimnames(bigarray) <- mydimnames
  save(bigarray, file = paste0('data/', "allsummaries.Rdata"))

  submethods <- mymethods
  for (num in c(20, 40, 65)) {
    load(file=paste0(outpath, "cparrays.RD.", num, ".",200,".Rdata"))
    bigarray[,submethods,,,paste(num), "RD"] <- arrays$summaries[,submethods,,,paste(num),]
    load(file=paste0(outpath, "cparrays.RR.", num, ".",200,".Rdata"))
    bigarray[,submethods,,,paste(num), "RR"] <- arrays$summaries[,submethods,,,paste(num),]
  }

  ### OPTIONAL: re-run calculations for large sample size
  # For RD
  system.time(cp205RD25 <- onecpfun(0.4, 0.1, n=205, contrast = "RD", alph=0.05, phis=0.25, methods=RDpairteam))[[3]]/60
  system.time(cp205RD75 <- onecpfun(0.3, 0.2, n=205, contrast = "RD", alph=0.05, phis=0.75, methods=RDpairteam))[[3]]/60
  system.time(cp205RD2599 <- onecpfun(0.4, 0.1, n=205, contrast = "RD", alph=0.01, phis=0.25, methods=RDpairteam))[[3]]/60
  system.time(cp205RD7599 <- onecpfun(0.3, 0.2, n=205, contrast = "RD", alph=0.01, phis=0.75, methods=RDpairteam))[[3]]/60
  system.time(cp205RD10 <- onecpfun(0.4, 0.1, n=205, contrast = "RD", alph=0.05, phis=0.10, methods=RDpairteam))[[3]]/60
  system.time(cp205RD50 <- onecpfun(0.4, 0.2, n=205, contrast = "RD", alph=0.05, phis=0.50, methods=RDpairteam))[[3]]/60
  system.time(cp205RD1099 <- onecpfun(0.4, 0.1, n=205, contrast = "RD", alph=0.01, phis=0.10, methods=RDpairteam))[[3]]/60
  system.time(cp205RD5099 <- onecpfun(0.4, 0.2, n=205, contrast = "RD", alph=0.01, phis=0.50, methods=RDpairteam))[[3]]/60
  # For RR
  system.time(cp205RR25 <- onecpfun(0.4, 0.1, n=205, contrast = "RR", alph=0.05, phis=0.25, methods=RRpairteam))[[3]]/60
  system.time(cp205RR75 <- onecpfun(0.3, 0.2, n=205, contrast = "RR", alph=0.05, phis=0.75, methods=RRpairteam))[[3]]/60
  system.time(cp205RR2599 <- onecpfun(0.4, 0.1, n=205, contrast = "RR", alph=0.01, phis=0.25, methods=RRpairteam))[[3]]/60
  system.time(cp205RR7599 <- onecpfun(0.3, 0.2, n=205, contrast = "RR", alph=0.01, phis=0.75, methods=RRpairteam))[[3]]/60
  system.time(cp205RR10 <- onecpfun(0.4, 0.1, n=205, contrast = "RR", alph=0.05, phis=0.10, methods=RDpairteam))[[3]]/60
  system.time(cp205RR50 <- onecpfun(0.4, 0.2, n=205, contrast = "RR", alph=0.05, phis=0.50, methods=RDpairteam))[[3]]/60
  system.time(cp205RR1099 <- onecpfun(0.4, 0.1, n=205, contrast = "RR", alph=0.01, phis=0.10, methods=RDpairteam))[[3]]/60
  system.time(cp205RR5099 <- onecpfun(0.4, 0.2, n=205, contrast = "RR", alph=0.01, phis=0.50, methods=RDpairteam))[[3]]/60
  # Combine results into a data object
  bignsummary <- array(NA, dim=c(6, 6, 4, 2, 2, 1))
  dimnames(bignsummary) <-
    c(dimnames(cp205RD10)[c(1,2)],
      list(c("0.4|0.1|0.1","0.4|0.1|0.25","0.4|0.2|0.5","0.3|0.2|0.75"),
           paste(c(0.05, 0.01)),
           c("RD", "RR"),
           "205"
      ))
  for (i in c(0.05, 0.01)) {
    for (k in contrasts) {
      bignsummary[,,"0.4|0.1|0.1", paste(i), paste(k), "205"] <-
        eval(parse(text=paste0("cp205", paste(k), "10", ifelse(i == 0.05, "", "99"))))
      bignsummary[,,"0.4|0.1|0.25", paste(i), paste(k), "205"] <-
        eval(parse(text=paste0("cp205", paste(k), "25", ifelse(i == 0.05, "", "99"))))
      bignsummary[,,"0.4|0.2|0.5", paste(i), paste(k), "205"] <-
        eval(parse(text=paste0("cp205", paste(k), "50", ifelse(i == 0.05, "", "99"))))
      bignsummary[,,"0.3|0.2|0.75", paste(i), paste(k), "205"] <-
        eval(parse(text=paste0("cp205", paste(k), "75", ifelse(i == 0.05, "", "99"))))
    }
  }
  save(bignsummary, file = paste(outpath, "bignsummary.Rdata"))


  ### FIGURE 1: CP, MACP, location index and DNCP for selected methods for RD, with N = 40, \alpha=0.05 and \phi=0.25
  load(file = paste0('data/', "cparrays.RD.", 40, ".",200,".Rdata"))
  plotpanel(plotdata = arrays, alpha = 0.05, par3 = 0.25,
            sel = c("SCAS-bc", "SCAS", "AS", "MOVER-NJ", "MOVER-W", "BP"),
            plotlab = "RDpair", fmt="tiff")

  ### FIGURE 2: CP, MACP, location index and DNCP for selected methods for RR, with N = 40, \alpha=0.05 and \phi=0.25
  load(file = paste0('data/', "cparrays.RR.", 40, ".",200,".Rdata"))
  plotpanel(plotdata = arrays, alpha = 0.05, par3 = 0.25,
            sel = c("SCAS-bc", "SCAS", "AS", "MOVER-NJ", "MOVER-W", "BP"),
            plotlab = "RRpair", fmt="tiff")

  ### FIGURE 3: Type I error for McNemar test and 'N-1' test
  # 2-D Type I error plot
  load(file=paste0(outpath, "cparrays.RD.", 65, ".",200,".Rdata"))
  p2 <- p1 <- seq(0, 1, length.out=201)
  cp1 <- onecpfun(
    p1 = p1,
    p2 = p2,
    ciarrays = arrays,
    alph = 0.05,
    phis = 0.25
  )
  res.factor <- 3
  tiff(file = paste0(outpath,"_tiff/","Laud_Fig3.tiff"),
       width=300*res.factor,
       height=300*res.factor,
       type="quartz"
  )
  par(pty='s')
  par(cex.main = res.factor*0.8*1, cex.axis=res.factor*0.8*1)
  par(mar = res.factor*(c(2,3,1,0.5)+0.1))
  plot(p2,
       1 - cp1[,"SCAS-bc","cp"],
       type = "l",
       ylim = c(0, 0.06),
       xlab = '',
       ylab = '',
       main = paste0("N = 65, \u03D5 = 0.25"),
       xaxt='n',
       yaxt='n',
       cex.lab = res.factor
  )
  axis(side = 2, las = 2)
  axis(side = 1, las = 1, )
  mtext(side = 1,
        text = bquote(paste(italic(p)[1]," = ",italic(p)[2])),
        cex = res.factor*1,
        line = 1.5*res.factor)
  mtext(side = 2,
        text = "Type I error rate",
        cex = res.factor*1,
        line = 2*res.factor)
  abline(h=0.05, lty=3, lwd=res.factor)
  lines(p2, 1 - cp1[,"SCAS-bc","cp"], lty=1, lwd = 2*res.factor)
  lines(p2, 1 - cp1[,"AS","cp"], lty=2, lwd = 2*res.factor)
  dev.off()

  ### TABLE 2.1 & TABLE 3.1: Summary of % proximate for selected methods for RD & RR
  load(file = paste0('data/', "allsummaries.Rdata"))
  mysummaries <-
    bigarray[, , c('95', '90'), c("meanCP", "pctCons", "pctnear", "pctgoodloc", "meanlocindex"),,]
  # Round to 0dps
  mysummaries[, , , c("pctCons", "pctnear", "pctgoodloc"),,] <-
    round(as.numeric(mysummaries[, , ,c("pctCons", "pctnear", "pctgoodloc"),,]), 0)
  # Add brackets to anticonservative methods
  anticons <- (as.numeric(mysummaries[,,,"pctCons",,]) < 50)
  mysummaries[,,,"pctnear",,] <- paste0(ifelse(anticons,"["," "),
                                        mysummaries[,,,"pctnear",,],
                                        ifelse(anticons,"]"," "))

  mytable1 <-
    ftable((mysummaries[, , , c("pctnear"),,]), col.vars = c(3,1), row.vars = c(5,4,2))
  write.ftable(mytable1, sep=',', quote=TRUE, file = paste0(outpath, "tablecp2065.csv"))


  ### TABLE 2.2 & TABLE 3.2: Summary of % central for selected methods for RD & RR
  mytable2 <-
    ftable((mysummaries[, , , c("pctgoodloc"),,]), col.vars = c(3,1), row.vars = c(5,4,2))
  write.ftable(mytable2, sep=',', quote=TRUE, file = paste0(outpath, "tableloc2065.csv"))


  ### Table 4: DNCP (One-sided type I error) for selected PSPs with larger sample size: N=205. Target DNCP=\ \alpha/2
  load(file=paste0('data/',"bignsummary.Rdata"))
  mytable2 <-
    ftable(bignsummary[,"dncp",,,,], col.vars = c(3,2), row.vars = c(4,1))
  write.ftable(round(mytable2, 4), sep=',', quote=TRUE, file = paste0(outpath, "bigntable205.csv"))


  ### Table 5.1: Example confidence intervals for RD, with (a, b, c, d) = (1, 1, 7, 12)
  x <- c(1, 1, 7, 12)
  egCI <- allpairci(x = x, contrast = "RD",
                    methods <- c("SCAS", "SCAS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP",
                                 "SCAS-cc125", "SCAS-cc5", "MOVER-NJcc125", "MOVER-NJcc5"),
                    alpha=0.05)
  dimnames(egCI)[[1]] <- ""
  ftable(round(egCI, 3), row.vars = 3)

  ### Table 5.2: Example confidence intervals for RR, with (a, b, c, d) = (1, 1, 7, 12)
  x <- c(1, 1, 7, 12)
  egCI <- allpairci(x = x, contrast = "RR",
                    methods <- c("SCAS", "SCAS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP",
                                 "SCAS-cc125", "SCAS-cc5", "MOVER-NJcc125", "MOVER-NJcc5"),
                    alpha=0.05)
  dimnames(egCI)[[1]] <- ""
  ftable(round(egCI, 3), row.vars = 3)


  ### p.5 footnote: MOVER-NJ intervals discrepancy vs M-L Tang et al
  xs <- rbind(c(43, 0, 1, 0),
             c(8, 3, 1, 2),
             c(4, 9, 3, 16))
  egCI <- allpairci(x = xs, contrast = "RD",
            methods <- c("AS", "MOVER-NW", "MOVER-NJ"),
            alpha=0.05)
  ftable(round(egCI, 4), row.vars = 1, col.vars=c(3, 2))


  ### p.12 evaluation of type I error rate (TIER) for test for association
  load(file=paste0('data/',"tiers.Rdata"))



  tier <- function(params) {
    n <- c(params[1])
    psi <- params[2]
    p1 <- p2 <- params[3]
    g <- expand.grid(x11 = 0:n, x12 = 0:n, x21 = 0:n)
    # reduce to possible combinations of a,b,c for paired data.
    g <- g[(g$x12 <= n - g$x11) &
             (g$x21 <= n - g$x11 - g$x12), ]
    xs <- data.matrix(cbind(g, x22 = n - rowSums(g)))
    prob <- pdfpair(p1 = p1,
                    p2 = p2,
                    psi = psi,
                    x = xs)
    xsub <- xs[prob > 1E-8, ]
    if (dim(xsub)[1] > 0) {
      pvals <- sapply(1:dim(xsub)[[1]], function(i)
        pchisq(scorepair(theta = 0,
                         x = xsub[i,],
                         contrast = "RD",
                         cc = FALSE,
                         skew = TRUE,
                         bcf = TRUE)$score^2, df=1, lower.tail=F) # 2-sided p-value
      )
      tier <- (pvals < 0.05) %*% prob[prob > 1E-8]
      pvals2 <- sapply(1:dim(xsub)[[1]], function(i)
        pchisq(scorepair(theta = 0,
                         x = xsub[i,],
                         contrast = "RD",
                         cc = FALSE,
                         skew = FALSE,
                         bcf = FALSE)$score^2, df=1, lower.tail=F)
      )
      tier2 <- (pvals2 < 0.05) %*% prob[prob > 1E-8]
    } else {
      tier <- 0
      tier2 <- 0
    }
    c(nminus1 = tier, mcnemar = tier2)
  }

  # Parameter scenarios matching Fagerland 2013
  myparams <- expand.grid(p1 = seq(0, 1, 0.01), psi = c(1, 2, 3, 5, 10), n = seq(10, 100, 5))
  # Runtime: 3.3 hours
  system.time(tiers <- sapply(1:dim(myparams)[1], function(i) tier(rev(unlist(myparams[i,])))))[[3]]/60
  # Alternative scenarios with larger correlations
#  myparams <- expand.grid(p1 = seq(0, 1, 0.02), phi = seq(0.25, 0.75, 0.05), n = seq(10, 100, 5))
#  system.time(tiers <- sapply(1:dim(myparams)[1], function(i) tier(rev(unlist(myparams[i,])))))[[3]]/60

  mytiers <- t(tiers)
  dimnames(mytiers)[[2]] <- c("SCAStier", "AStier")
  mytiers <- cbind(myparams, mytiers)
  save(mytiers, file = paste0(outpath, "tiers.Rdata"))

  apply(mytiers[,4:5], 2, mean)
  apply(mytiers[,4:5], 2, max)
  apply(mytiers[,4:5], 2, function(x) mean(x > 0.05))
  apply(mytiers[,4:5], 2, function(x) mean(x < 0.03))
  dim(mytiers)


}
