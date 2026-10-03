set.seed(2012) #ensure we use the same jitters for each run

if (FALSE) {
  # The following code reproduces the results in the manuscript
  # 'Equal-tailed confidence intervals for paired binomial proportions' by Peter J. Laud
  #
  # Open the cpplot.Rproject file, then:
  # > devtools::load.all(".")
  # then run selected code from within 'create_outputs.R'

  ### NOTE methods have been re-labelled in the manuscript:
  ### SCAS --> 'SCASu' (method without the 'N-1' variance bias correction)
  ### SCAS-bc --> 'SCAS' (including the 'N-1' bias correction)

  # Set path for output files as required by user
  outpath <- '/myoutputpath/'
  outpath <- paste0(root, "Main/Courses_papers/skewscore/paired/") # Remove for final upload
  outpath <- "D:/Pete/Documents/Research/paired/" # Remove for final upload
  outpath1 <- "D:/Pete/Documents/Research/paired/" # Remove for final upload
  outpath <- "D:/Pete/Documents/Research/paired/negcorr/" # Remove for final upload

  #  outpath <- 'data/'

  #############################################################################
  ### OPTIONAL (summary data from these runs is provided in the repository)
  ### Run the CP calculation function for N=20, N=40 and N=65
  ### WARNING: for N=40 and 65, these take several hours to run!
  #############################################################################
  RDpairteam <- c("SCAS-bc", "AS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP", "AM", "Wald")
  RRpairteam <- c("SCAS-bc", "AS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP", "BP-J", "Wald")
  # ^^ Example of how you might specify a subset of methods used for larger N
  # for reduced runtimes using methods= argument in cpfun() below
  alphas <- c(0.1, 0.05, 0.01)
  phis <- c(-0.1, 0.1, 0.25, 0.5, 0.75)

  system.time(mycis <- cifun(n=20, contrast="RD", alph = alphas))[[3]]/60
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
  system.time(mycis <- cifun(n=20, contrast="RR", alph = alphas))[[3]]/60
# load(file=paste0(outpath1, "cis.RR.20.Rdata")); mycis <- ciarrays
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
  system.time(mycis <- cifun(n=40, contrast="RD", alph = alphas))[[3]]/60
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
  system.time(mycis <- cifun(n=40, contrast="RR", alph = alphas))[[3]]/60
  # load(file=paste0(outpath1, "cis.RR.40.Rdata")); mycis <- ciarrays
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
  system.time(mycis <- cifun(n=65, contrast="RD", alph = alphas))[[3]]/60
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
  system.time(mycis <- cifun(n=65, contrast="RR", alph = alphas))[[3]]/60
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60

  # Evaluation for conditional OR
  ORpairteam <- c("SCASp", "SCASpu", "mid-p", "Jeffreys", "Wilson", "Blaker", "Wald", "Laplace")
  #                "SCASp-c125", "SCASp-c25", "SCASp-c5", "midp-c25", "Jeffreys-c125", "Jeffreys-c25", "C-P") # Subset of OR methods used for larger N
  system.time(mycis <- cifun(n=20, contrast="OR", alph = alphas))[[3]]/60
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
  system.time(mycis <- cifun(n=40, contrast="OR", alph = alphas))[[3]]/60
  # load(file=paste0(outpath1, "cis.OR.40.Rdata")); mycis <- ciarrays
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
  system.time(mycis <- cifun(n=65, contrast="OR", alph = alphas))[[3]]/60
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
  # larger sample evaluation for OR with coarser grid of PSPs
  # - Note: OR methods with closed form expressions are quicker to calculate,
  #   but the coverage probability calculations take longer because there is an extra step
  #   to get p12, p21 from p1, p2 and phi
  alphas <- 0.05
#  system.time(mycis <- cifun(n=105, contrast="OR", alph = alphas, methods = ORpairteam))[[3]]/60
#  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=100, jitt=T, smooth=T, phis=phis))[[3]]/60

  RDmeth <- c("SCAS-bc", "SCAS", "AS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP", "Wald", "Wald-cc",
              "SCAS-c5", "SCAS-c25", "SCAS-c125", "MOVER-c5", "MOVER-c25", "MOVER-c125")
  RRmeth <- c("SCAS-bc", "SCAS", "AS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP", "BP-J", "Wald",
              "SCAS-c5", "SCAS-c25", "SCAS-c125", "MOVER-c5", "MOVER-c25", "MOVER-c125", "BP-cc")
  ORmeth <- c("SCASp", "SCASpu", "mid-p", "Jeffreys", "Wilson", "Blaker", "Wald", "Laplace",
              "SCASp-c5", "SCASp-c25", "SCASp-c125", "C-P", "Jeffreys-c25", "Jeffreys-c125")


  # Limited versions for GitHub due to file size limit
  load(file=paste(outpath, "cparrays.RD.", 40, ".",200,".Rdata",sep=""))
  mycis <- arrays
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, alph=0.05, phis=phis,
                                          methods = RDmeth[1:6], outdir="data/"))[[3]]/60
  load(file=paste(outpath, "cparrays.RR.", 40, ".",200,".Rdata",sep=""))
  mycis <- arrays
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, alph=0.05, phis=phis,
                                          methods = RRmeth[1:7], outdir="data/"))[[3]]/60
  load(file=paste(outpath, "cparrays.OR.", 40, ".",200,".Rdata",sep=""))
  mycis <- arrays
  Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, alph=0.05, phis=phis,
                                          methods = ORmeth[1:7], outdir="data/"))[[3]]/60

  # Subset arrays to selected N and methods for smaller file size to upload to GitHub
  load(file=paste(outpath, "cparrays.RD.", 40, ".",200,".Rdata",sep=""))
  arrays$cis <- arrays$cis[,, RDmeth[1:6], "95",,, drop=FALSE]
  arrays$mastercp <- arrays$mastercp[,,, RDmeth[1:6], "95",,,, drop=FALSE]
  save(arrays, file = paste("data/cparrays.RD.40.200.Rdata", sep = ""))

  load(file=paste(outpath, "cparrays.RR.", 40, ".",200,".Rdata",sep=""))
  arrays$cis <- arrays$cis[,, RRmeth[1:7], "95",,, drop=FALSE]
  arrays$mastercp <- arrays$mastercp[,,, RDmeth[1:7], "95",,,, drop=FALSE]
  save(arrays, file = paste("data/cparrays.RR.40.200.Rdata", sep = ""))

  load(file=paste(outpath, "cparrays.OR.", 40, ".",200,".Rdata",sep=""))
  arrays$cis <- arrays$cis[,, ORmeth[1:5], "95",,, drop=FALSE]
  arrays$mastercp <- arrays$mastercp[,,, ORmeth[1:5], "95",,,, drop=FALSE]
  save(arrays, file = paste("data/cparrays.OR.40.200.Rdata", sep = ""))


  #############################################################################
  # OPTIONAL:
  # Combine output arrays for summarising across different Ns
  # Note these are not all included in the GitHub repository due to size
  # So the following requires all the cpfun calls above to be run
  #############################################################################
  mynums <- c(20, 40, 65)
  phis <- c(0.1, 0.25, 0.5, 0.75)
  mymethods <- c("SCAS-bc", "SCAS", "AS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W",
                 "AM", "BP", "BP-J", "BP-cc", "Wald", "Wald-cc",
                 "SCAS-c5", "SCAS-c25", "SCAS-c125", "MOVER-c5", "MOVER-c25", "MOVER-c125",
                 "SCASp", "SCASpu", "Jeffreys", "mid-p", "Wilson", "Blaker", "Laplace",
                 "SCASp-c5", "SCASp-c25", "SCASp-c125", "C-P", "Jeffreys-c25", "Jeffreys-c125")
  nmeth <- length(mymethods)
  load(file=paste0(outpath, "cparrays.RD.", 40, ".",200,".Rdata"))
  mydims <- dim(arrays$summaries)
  mydims[5] <- length(mynums)
  mydims[6] <- 3
  mydims[2] <- nmeth
  mydimnames <- dimnames(arrays$summaries)
  mydimnames[[5]] <- paste(mynums)
  mydimnames[[6]] <- c("RD", "RR", "OR")
  mydimnames[[2]] <- mymethods

  bigarray <- array(NA, dim = mydims)
  dimnames(bigarray) <- mydimnames

  for (num in mynums) {
    load(file=paste0(outpath, "cparrays.RD.", num, ".",200,".Rdata"))
    bigarray[,RDmeth,,,paste(num), "RD"] <- arrays$summaries[,RDmeth,,,paste(num),]
  }
  for (num in mynums) {
    load(file=paste0(outpath, "cparrays.RR.", num, ".",200,".Rdata"))
    bigarray[,RRmeth,,,paste(num), "RR"] <- arrays$summaries[,RRmeth,,,paste(num),]
  }
  for (num in mynums) {
    load(file=paste0(outpath, "cparrays.OR.", num, ".",200,".Rdata"))
    bigarray[,ORmeth,,,paste(num), "OR"] <- arrays$summaries[paste(phis), ORmeth,,,paste(num),]
  }


  save(bigarray, file = paste0(outpath, "allsummaries.Rdata"))


  #############################################################################
  ### OPTIONAL: re-run calculations for large sample size (Table 5) (takes several hours)
  #############################################################################
  RDpairteam <- c("SCAS-bc", "SCAS", "AS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP", "AM", "Wald")
#                  "Wald-cc", "SCAS-c125") # Waldcc included to match vector length for RR
  RRpairteam <- c("SCAS-bc", "SCAS", "AS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP", "BP-J", "Wald")
  # For RD
  system.time(cp205RD10 <- onecpfun(p1=0.1, p2=0.4, n=205, contrast = "RD", alph=0.05, phis=0.10, methods=RDpairteam))[[3]]/60
  system.time(cp205RD25 <- onecpfun(0.1, 0.4, n=205, contrast = "RD", alph=0.05, phis=0.25, methods=RDpairteam))[[3]]/60
  system.time(cp205RD50 <- onecpfun(0.2, 0.4, n=205, contrast = "RD", alph=0.05, phis=0.50, methods=RDpairteam))[[3]]/60
  system.time(cp205RD75 <- onecpfun(0.3, 0.4, n=205, contrast = "RD", alph=0.05, phis=0.75, methods=RDpairteam))[[3]]/60
  system.time(cp205RD95 <- onecpfun(0.39, 0.4, n=205, contrast = "RD", alph=0.05, phis=0.95, methods=RDpairteam))[[3]]/60
  system.time(cp205RD1099 <- onecpfun(0.1, 0.4, n=205, contrast = "RD", alph=0.01, phis=0.10, methods=RDpairteam))[[3]]/60
  system.time(cp205RD2599 <- onecpfun(0.1, 0.4, n=205, contrast = "RD", alph=0.01, phis=0.25, methods=RDpairteam))[[3]]/60
  system.time(cp205RD5099 <- onecpfun(0.2, 0.4, n=205, contrast = "RD", alph=0.01, phis=0.50, methods=RDpairteam))[[3]]/60
  system.time(cp205RD7599 <- onecpfun(0.3, 0.4, n=205, contrast = "RD", alph=0.01, phis=0.75, methods=RDpairteam))[[3]]/60
  system.time(cp205RD9599 <- onecpfun(0.39, 0.4, n=205, contrast = "RD", alph=0.01, phis=0.95, methods=RDpairteam))[[3]]/60
  # For RR
  system.time(cp205RR10 <- onecpfun(0.1, 0.4, n=205, contrast = "RR", alph=0.05, phis=0.10, methods=RRpairteam))[[3]]/60
  system.time(cp205RR25 <- onecpfun(0.1, 0.4, n=205, contrast = "RR", alph=0.05, phis=0.25, methods=RRpairteam))[[3]]/60
  system.time(cp205RR50 <- onecpfun(0.2, 0.4, n=205, contrast = "RR", alph=0.05, phis=0.50, methods=RRpairteam))[[3]]/60
  system.time(cp205RR75 <- onecpfun(0.3, 0.4, n=205, contrast = "RR", alph=0.05, phis=0.75, methods=RRpairteam))[[3]]/60
  system.time(cp205RR95 <- onecpfun(0.39, 0.4, n=205, contrast = "RR", alph=0.05, phis=0.95, methods=RRpairteam))[[3]]/60
  system.time(cp205RR1099 <- onecpfun(0.1, 0.4, n=205, contrast = "RR", alph=0.01, phis=0.10, methods=RRpairteam))[[3]]/60
  system.time(cp205RR2599 <- onecpfun(0.1, 0.4, n=205, contrast = "RR", alph=0.01, phis=0.25, methods=RRpairteam))[[3]]/60
  system.time(cp205RR5099 <- onecpfun(0.2, 0.4, n=205, contrast = "RR", alph=0.01, phis=0.50, methods=RRpairteam))[[3]]/60
  system.time(cp205RR7599 <- onecpfun(0.3, 0.4, n=205, contrast = "RR", alph=0.01, phis=0.75, methods=RRpairteam))[[3]]/60
  system.time(cp205RR9599 <- onecpfun(0.39, 0.4, n=205, contrast = "RR", alph=0.01, phis=0.95, methods=RRpairteam))[[3]]/60
  # Combine results into a data object
  bignsummary <- array(NA, dim=c(10, 6, 5, 2, 2, 1))
  dimnames(bignsummary) <-
    list(RRpairteam, unlist(dimnames(cp205RD75)[2]),
      c("0.1|0.4|0.1","0.1|0.4|0.25","0.2|0.4|0.5","0.3|0.4|0.75","0.39|0.4|0.95"),
           paste(c(0.05, 0.01)),
           c("RD", "RR"),
           "205"
      )
  for (i in c(0.05, 0.01)) {
    for (k in c("RD", "RR")) {
      teamsel <- eval(parse(text=paste0(k, "pairteam")))
      bignsummary[teamsel,,"0.1|0.4|0.1", paste(i), paste(k), "205"] <-
        eval(parse(text=paste0("cp205", paste(k), "10", ifelse(i == 0.05, "", "99"))))[teamsel,,]
      bignsummary[teamsel,,"0.1|0.4|0.25", paste(i), paste(k), "205"] <-
        eval(parse(text=paste0("cp205", paste(k), "25", ifelse(i == 0.05, "", "99"))))[teamsel,,]
      bignsummary[teamsel,,"0.2|0.4|0.5", paste(i), paste(k), "205"] <-
        eval(parse(text=paste0("cp205", paste(k), "50", ifelse(i == 0.05, "", "99"))))[teamsel,,]
      bignsummary[teamsel,,"0.3|0.4|0.75", paste(i), paste(k), "205"] <-
        eval(parse(text=paste0("cp205", paste(k), "75", ifelse(i == 0.05, "", "99"))))[teamsel,,]
      bignsummary[teamsel,,"0.39|0.4|0.75", paste(i), paste(k), "205"] <-
        eval(parse(text=paste0("cp205", paste(k), "75", ifelse(i == 0.05, "", "99"))))[teamsel,,]
    }
  }
  save(bignsummary, file = paste(outpath, "bignsummary.Rdata"))

  # For OR
  ORpairteam <- c("SCASp", "SCASpu", "mid-p", "Jeffreys", "Wilson", "Blaker", "Wald", "Laplace")
  system.time(cp205OR10 <- onecpfun(0.1, 0.4, n=205, contrast = "OR", alph=0.05, phis=0.10, methods=ORpairteam))[[3]]/60
  system.time(cp205OR25 <- onecpfun(0.1, 0.4, n=205, contrast = "OR", alph=0.05, phis=0.25, methods=ORpairteam))[[3]]/60
  system.time(cp205OR50 <- onecpfun(0.2, 0.4, n=205, contrast = "OR", alph=0.05, phis=0.50, methods=ORpairteam))[[3]]/60
  system.time(cp205OR75 <- onecpfun(0.3, 0.4, n=205, contrast = "OR", alph=0.05, phis=0.75, methods=ORpairteam))[[3]]/60
  system.time(cp205OR95 <- onecpfun(0.39, 0.4, n=205, contrast = "OR", alph=0.05, phis=0.95, methods=ORpairteam))[[3]]/60
  system.time(cp205OR1099 <- onecpfun(0.1, 0.4, n=205, contrast = "OR", alph=0.01, phis=0.10, methods=ORpairteam))[[3]]/60
  system.time(cp205OR2599 <- onecpfun(0.1, 0.4, n=205, contrast = "OR", alph=0.01, phis=0.25, methods=ORpairteam))[[3]]/60
  system.time(cp205OR5099 <- onecpfun(0.2, 0.4, n=205, contrast = "OR", alph=0.01, phis=0.50, methods=ORpairteam))[[3]]/60
  system.time(cp205OR7599 <- onecpfun(0.3, 0.4, n=205, contrast = "OR", alph=0.01, phis=0.75, methods=ORpairteam))[[3]]/60
  system.time(cp205OR9599 <- onecpfun(0.39, 0.4, n=205, contrast = "OR", alph=0.01, phis=0.95, methods=ORpairteam))[[3]]/60
  bignsummaryOR <- array(NA, dim=c(8, 6, 5, 2, 1, 1))
  dimnames(bignsummaryOR) <-
    c(dimnames(cp205OR10)[c(1,2)],
      list(c("0.1|0.4|0.1","0.1|0.4|0.25","0.2|0.4|0.5","0.3|0.4|0.75","0.39|0.4|0.95"),
           paste(c(0.05, 0.01)),
           c("OR"),
           "205"
      ))
  for (i in c(0.05, 0.01)) {
    for (k in "OR") {
      bignsummaryOR[,,"0.4|0.1|0.1", paste(i), paste(k), "205"] <-
        eval(parse(text=paste0("cp205", paste(k), "10", ifelse(i == 0.05, "", "99"))))
      bignsummaryOR[,,"0.4|0.1|0.25", paste(i), paste(k), "205"] <-
        eval(parse(text=paste0("cp205", paste(k), "25", ifelse(i == 0.05, "", "99"))))
      bignsummaryOR[,,"0.4|0.2|0.5", paste(i), paste(k), "205"] <-
        eval(parse(text=paste0("cp205", paste(k), "50", ifelse(i == 0.05, "", "99"))))
      bignsummaryOR[,,"0.3|0.2|0.75", paste(i), paste(k), "205"] <-
        eval(parse(text=paste0("cp205", paste(k), "75", ifelse(i == 0.05, "", "99"))))
    }
  }
  save(bignsummaryOR, file = paste0(outpath, "bignsummaryOR.Rdata"))


  #############################################################################
  ### FIGURE 1: 2-D plots of CP and LNCP for selected methods. See cpslice.R
  #############################################################################



  #############################################################################
  ### FIGURE 2: CP, MACP, LNCP, location index and DNCP for selected methods for RD, with N = 40, \alpha=0.05 and \phi=0.25
  #############################################################################
  load(file = paste0(outpath, "cparrays.RD.", 40, ".",200,".Rdata"))
  plotpanel(plotdata = arrays, alpha = 0.05, par3 = 0.25,
            sel = c("SCAS-bc", "AS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP", "Wald"),
            plotlab = "RDpairW", fmt="tiff", res.factor = 6, CIlen = TRUE)

  #############################################################################
  ### FIGURE 3: CP, MACP, location index and DNCP for selected methods for RR, with N = 40, \alpha=0.05 and \phi=0.25
  #############################################################################
  load(file = paste0(outpath, "cparrays.RR.", 40, ".", 200, ".Rdata"))
  plotpanel(plotdata = arrays, alpha = 0.05, par3 = 0.25,
            sel = c("SCAS-bc", "AS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP", "Wald"),
            plotlab = "RRpairW", fmt="tiff", res.factor = 6, CIlen = TRUE)

  #############################################################################
  ### FIGURE 4: CP, MACP, location index and DNCP for selected methods for OR, with N = 40, \alpha=0.05 and \phi=0.25
  #############################################################################
  load(file = paste0(outpath, "cparrays.OR.", 40, ".",100,".Rdata"))
#  dimnames(arrays$summaries)
  plotpanel(plotdata = arrays, alpha = 0.05, par3 = 0.25,
            sel = c("SCASp", "SCASpu", "mid-p", "Jeffreys", "Wilson", "Wald", "Laplace"),
            plotlab = "ORpairW", fmt="tiff", res.factor = 6, CIlen = TRUE)

  #############################################################################
  ### FIGURE 5: Type I error for McNemar tests, mid-p test and 'N-1' AS test
  ### Uses TIERs dataset - see further down
  #############################################################################

  # Load the saved dataset
  load(file=paste0(outpath,"newtiers1.Rdata"))

  # Create a plot of TIERs
  mytiers <- tiers1
  res.factor <- 3
  grid.factor <- 2
  tiff(file = paste0(outpath,"_tiff/","Laud_Fig5new.tiff"),
       width = 300*grid.factor*res.factor,
       height = 600*res.factor,
       type = "windows"
       #       type="quartz"
  )
  #  par(pty='s')
  par(mfrow = c(2, 2))
  par(cex.main = grid.factor*res.factor*0.8*1, cex.axis=grid.factor*res.factor*0.5*1)
  #  par(mar = res.factor*(c(2,3,3,0.5)+0.1))
  methods <- c("nminus1", "midp", "mcnemar", "mcnemarcc")
  labels <- c("\'N - 1\' AS", "mid-p", "McNemar asymptotic", "McNemar asymptotic (cc)")
  for (i in 1:4) {
    par(mar = grid.factor*res.factor*(c(2,3,3,0.5)+0.1))
    plot(mytiers$p1,
         eval(parse(text=paste0("mytiers$", methods[i]))),
         type = "n",
         ylim = c(0, 0.06),
         xlab = '',
         ylab = '',
         main = labels[i],
         xaxt='n',
         yaxt='n',
         cex.lab = res.factor
    )
    axis(side = 2, las = 2)
    axis(side = 1, las = 1, padj=1)
    mtext(side = 1,
          text = bquote(paste(italic(p)[1]," = ",italic(p)[2])),
          cex = res.factor*1,
          line = 1.5*1.5*res.factor)
    mtext(side = 2,
          text = "Type I error rate",
          cex = res.factor*1,
          line = 1.5*2*res.factor)
    abline(h=0.05, lty=3, lwd=res.factor)

    #    for (ps in c(1, 2, 3, 5, 10)) {
    for (ps in unique(mytiers[,3])) {
      for (n in nseq2[1:4]) {
        #        tiersub <- mytiers[mytiers$psi == ps & mytiers$n == n, ]
        tiersub <- mytiers[mytiers[,3] == ps & mytiers$n == n, ]
        lines(tiersub$p1,
              eval(parse(text=paste0("tiersub$", methods[i]))),
              lty = 2,
              lwd = 0.5*res.factor,
              col = "gray50")
      }
      for (n in nseq2[5:length(nseq2)]) {
        #        tiersub <- mytiers[mytiers$psi == ps & mytiers$n == n, ]
        tiersub <- mytiers[mytiers[,3] == ps & mytiers$n == n, ]
        lines(tiersub$p1,
              eval(parse(text=paste0("tiersub$", methods[i]))),
              lty = 1,
              lwd = 0.25*res.factor)
      }
    }

  }
  dev.off()





  #############################################################################
  ### FIGURE 6: CP, MACP, location index and DNCP for selected conservative methods for RD, with N = 40, \alpha=0.05 and \phi=0.25
  #############################################################################
  load(file = paste0(outpath, "cparrays.RD.", 40, ".",200,".Rdata"))
  plotpanel(plotdata = arrays, alpha = 0.05, par3 = 0.25,
            sel = c("SCAS-c125", "SCAS-c25", "SCAS-c5", "MOVER-c125", "MOVER-c25", "MOVER-c5", "Wald-cc"),
            plotlab = "RDcpair", fmt="tiff", res.factor = 6, CIlen = TRUE)


  #############################################################################
  ### FIGURE 7: CP, MACP, location index and DNCP for selected conservative methods for RR, with N = 40, \alpha=0.05 and \phi=0.25
  #############################################################################
  load(file = paste0(outpath, "cparrays.RR.", 40, ".", 200, ".Rdata"))
  plotpanel(plotdata = arrays, alpha = 0.05, par3 = 0.25,
            sel = c("SCAS-c125", "SCAS-c25", "SCAS-c5", "MOVER-c125", "MOVER-c25", "MOVER-c5", "BP-cc"),
            plotlab = "RRcpair", fmt="tiff", res.factor = 6, CIlen = TRUE)

  #############################################################################
  ### FIGURE 8: CP, MACP, location index and DNCP for selected methods for OR, with N = 40, \alpha=0.05 and \phi=0.25
  #############################################################################
  load(file = paste0(outpath, "cparrays.OR.", 40, ".",200,".Rdata"))
  plotpanel(plotdata = arrays, alpha = 0.05, par3 = 0.25,
            sel = c("SCASp-c125", "SCASp-c25", "SCASp-c5", "midp-c25", "Jeffreys-c125", "Jeffreys-c25", "C-P", "Blaker"),
            plotlab = "ORcpair", fmt="tiff", res.factor = 6, CIlen = TRUE)


  #############################################################################
  ### TABLE 2, 3 & 4: Summary of each metric for selected methods for RD & RR, & OR
  #############################################################################
  load(file = paste0(outpath, "allsummaries.Rdata"))
  mysummaries <-
    bigarray[, , c('95', '90'), c("meanCP", "pctCons", "pctnear", "pctnear.1side", "pctnear.DNCP",
                                  "pctgoodloc", "meanlocindex", "pctBad.DNCP", "minCP", "pctCons.both"),
             c("20", "40", "65"), ,drop=F]
  # Round to 0dps
  mysummaries[, , , c("pctCons", "pctnear", "pctgoodloc", "pctBad.DNCP", "pctnear.1side", "pctnear.DNCP"),,] <-
    round(as.numeric(mysummaries[, , ,c("pctCons", "pctnear", "pctgoodloc", "pctBad.DNCP", "pctnear.1side", "pctnear.DNCP"),,]), 0)
  # Add brackets to anticonservative methods
  anticons <- (as.numeric(mysummaries[,,,"pctCons",,]) < 50)
  mysummaries[,,,"pctnear",,] <- paste0(ifelse(anticons,"["," "),
                                        mysummaries[,,,"pctnear",,],
                                        ifelse(anticons,"]"," "))

  dimnames(mysummaries)
#  RDmeth <- c("SCAS-bc", "SCAS", "AS", "MOVER-NJ", "MOVER-W", "BP")
#  RRmeth <- c("SCAS-bc", "SCAS", "AS", "MOVER-NJ", "MOVER-W", "BP", "BP-J")
#  ORmeth <- c("SCASp", "SCASpu", "mid-p", "Jeffreys", "Wilson")

  # Summarise by metric for N=40 for new tables in resubmission
  mytable2a <- ftable((mysummaries[, RDmeth[c(1, 3:9)], ,
                                   #                      c("pctnear", "pctnear.1side", "pctnear.DNCP", "pctBad.DNCP", "meanlocindex"),
                                   c("pctnear", "pctnear.1side", "pctBad.DNCP", "meanlocindex"),
                                   "40",c("RD")]), col.vars = c(3,1), row.vars = c(4,2))
  mytable2b <- ftable((mysummaries[, RDmeth[c(1, 3:9)], ,
                                   #                      c("pctnear", "pctnear.1side", "pctnear.DNCP", "pctBad.DNCP", "meanlocindex"),
                                   c("pctnear", "pctnear.1side", "pctBad.DNCP", "meanlocindex"),
                                   "65",c("RD")]), col.vars = c(3,1), row.vars = c(4,2))
  write.ftable(mytable2a, sep=',', quote=TRUE, justify="none", file = paste0(outpath, "table2a.csv"))
  write.ftable(mytable2b, sep=',', quote=TRUE, justify="none", file = paste0(outpath, "table2b.csv"))

  ftable((mysummaries[, RRmeth[c(1, 3:10)], ,
#                      c("pctnear", "pctnear.1side", "pctnear.DNCP", "pctBad.DNCP", "meanlocindex"),
                      c("pctnear", "pctnear.1side", "pctBad.DNCP", "meanlocindex"),
                      "40", c("RR")]), col.vars = c(3,1), row.vars = c(4,2))
  ftable((mysummaries[, ORmeth[c(1:5, 7:8)], ,
#                      c("pctnear", "pctnear.1side", "pctnear.DNCP", "pctgoodloc", "pctBad.DNCP", "meanlocindex"),
                      c("pctnear", "pctnear.1side", "pctBad.DNCP", "meanlocindex"),
                      "40",c("OR")]), col.vars = c(3,1), row.vars = c(4,2))


  #############################################################################
  ### Table 5: DNCP (One-sided type I error) for selected PSPs with larger sample size: N=205. Target DNCP=\ \alpha/2
  #############################################################################
  load(file=paste0(outpath,"bignsummary.Rdata"))
  mytable5aa <- round(bignsummary[,"dncp",,,,], 4)
  limit <- mytable5aa
  limit[,,"0.01", ] <- 0.0055
  limit[,,"0.05", ] <- 0.0275
  # Add brackets to methods with unacceptable coverage
  anticons <- (mytable5aa > limit)
  anticons[is.na(anticons)] <- FALSE
  mytable5at <- apply(mytable5aa, 1:4,
                      function(x) {
                              formatC(round(x, 4), format='f', digits=4)
                      }
                      )
  mytable5at[anticons] <- paste0("[", mytable5at[anticons], "]")


  mytable5a <-
#    round(ftable(bignsummary[,"dncp",,,,], col.vars = c(3,2), row.vars = c(4,1)), 4)
    ftable(mytable5at, col.vars = c(3,2), row.vars = c(4,1))
  write.ftable(mytable5a, sep=',', quote=TRUE, justify="none", file = paste0(outpath, "bigntable205.csv"))


  load(file=paste0(outpath,"bignsummaryOR.Rdata"))
  mytable5ba <- round(bignsummaryOR[,"dncp",,,,], 4)
  dim(mytable5ba)
  limit <- mytable5ba
  dimnames(limit)
  limit[,,"0.01"] <- 0.0055
  limit[,,"0.05"] <- 0.0275
  # Add brackets to methods with unacceptable coverage
  anticons <- (mytable5ba > limit)
  anticons[is.na(anticons)] <- FALSE
  mytable5bt <- apply(mytable5ba, 1:3,
                      function(x) {
                        formatC(round(x, 4), format='f', digits=4)
                      }
  )
  mytable5bt[anticons] <- paste0("[", mytable5bt[anticons], "]")



  mytable5b <-
#    ftable(bignsummaryOR[,"dncp",,,,], col.vars = c(3,2), row.vars = c(1))
    ftable(mytable5bt, col.vars = c(3,2), row.vars = c(1))
  write.ftable(mytable5b, sep=',', quote=TRUE, justify="none", file = paste0(outpath, "bignORtable205.csv"))


  #############################################################################
  ### Table 6: Conservative coverage summary
  #############################################################################

  RDcpairteam <- c("SCAS-c5", "SCAS-c25", "SCAS-c125", "SCAS", "MOVER-c5", "MOVER-c25", "MOVER-c125", "MOVER-NJ", "Wald-cc") 	#Paired RD, cc
  RRcpairteam <- c("SCAS-c5", "SCAS-c25", "SCAS-c125", "SCAS", "MOVER-c5", "MOVER-c25", "MOVER-c125", "MOVER-NJ", "BP-cc") 	#Paired RR, cc
  ORcpairteam <- c("SCASp-c5", "SCASp-c25", "SCASp-c125", "SCASp", "C-P", "Blaker", "Jeffreys-c25", "Jeffreys-c125", "Jeffreys")

  # Overall minimum coverage for each continuity-adjusted method per contrast
  # (Including corresponding unadjusted method for reference)
  apply(bigarray[,RDcpairteam,c("95","90"),"minCP",,"RD"], 2:3, min)
  apply(bigarray[,RRcpairteam,c("95","90"),"minCP",,"RR"], 2:3, min)
  apply(bigarray[,ORcpairteam,c("95","90"),"minCP",,"OR"], 2:3, min)

  # Overall average of %PSP that are doubly conservative (showing 2dps for results close to 100)
  apply(bigarray[,RDcpairteam,c("95","90"),"pctCons.both",,"RD"], 2:3,  function(x) round(mean(as.numeric(x)), 2))
  apply(bigarray[,RDcpairteam,c("95","90"),"pctCons.both",,"RD"], 2:3,  function(x) round(mean(as.numeric(x)), 0))
  apply(bigarray[,RRcpairteam,c("95","90"),"pctCons.both",,"RR"], 2:3,  function(x) round(mean(as.numeric(x)), 2))
  apply(bigarray[,RRcpairteam,c("95","90"),"pctCons.both",,"RR"], 2:3,  function(x) round(mean(as.numeric(x)), 0))
  apply(bigarray[,ORcpairteam,c("95","90"),"pctCons.both",,"OR"], 2:3,  function(x) round(mean(as.numeric(x)), 2))
  apply(bigarray[,ORcpairteam,c("95","90"),"pctCons.both",,"OR"], 2:3,  function(x) round(mean(as.numeric(x)), 0))

  # Experimental: Overall average location index for N>=40
  apply(bigarray[,RDpairteam,c("95","90"),"meanlocindex",c(2:3),"RD"], 2:3,  function(x) round(mean(as.numeric(x)), 2))
  apply(bigarray[,RRpairteam,c("95","90"),"meanlocindex",2:3,"RR"], 2:3,  function(x) round(mean(as.numeric(x)), 2))
  apply(bigarray[,ORpairteam,c("95","90"),"meanlocindex",2:3,"OR"], 2:3,  function(x) round(mean(as.numeric(x)), 2))



  #############################################################################
  ### Table 7: Example confidence intervals with (a, b, c, d) = (1, 1, 7, 12)
  #############################################################################
  x <- c(1, 1, 7, 12)
  egCI <- allpairci(x = x, contrast = "RD",
                    methods = c("SCAS-bc", "AS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP",
                                 "SCAS-c125", "SCAS-c5", "AS-bc-c125", "AS-bc-c5", "MOVER-c125", "MOVER-c5"),
                    alpha=0.05)
  dimnames(egCI)[[1]] <- ""
  ftable(round(egCI, 3), row.vars = 3)
  # Calculate width from unrounded data, for consistency with Fagerland et al
  t(t((round(egCI[,2,], 3) - round(egCI[,1,],3))))

  egCI <- allpairci(x = x, contrast = "RR",
                    methods = c("SCAS-bc", "AS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP", "BP-J",
                                 "SCAS-c125", "SCAS-c5", "AS-bc-c125", "AS-bc-c5", "MOVER-c125", "MOVER-c5"),
                    alpha=0.05)
  dimnames(egCI)[[1]] <- ""
  ftable(round(egCI, 3), row.vars = 3)
  # Calculate log width from unrounded data, for consistency with Fagerland et al
  t(t((round(log(egCI[,2,]) - log(egCI[,1,]),2))))


  egCI <- allpairci(x = x, contrast = "OR",
                    methods = c("SCASp", "SCASpu", "Wilson", "Jeffreys", "mid-p", "Wald", "Laplace",
                                 "SCASp-c125", "SCASp-c5", "Jeffreys-c125", "C-P"),
                    alpha=0.05)
  dimnames(egCI)[[1]] <- ""
  ftable(round(egCI, 3), row.vars = 3)
  t(t((round(log(egCI[,2,]) - log(egCI[,1,]),2))))
  (round(log(egCI[,2,]) - log(egCI[,1,]),2))


  #############################################################################
  ### p.6 footnote: MOVER-NJ RD intervals discrepancy vs M-L Tang et al 2020a Table VI
  ### (cf methods AS=TANGO, MOVER-NW=NW and MOVER=NJ≠NJ)
  #############################################################################
  xs <- rbind(c(43, 0, 1, 0),
             c(8, 3, 1, 2),
             c(4, 9, 3, 16))
  egCI <- allpairci(x = xs, contrast = "RD",
            methods <- c("AS", "MOVER-NW", "MOVER-NJ"),
            alpha=0.05)
  ftable(round(egCI, 4), row.vars = 1, col.vars=c(3, 2))

  # And RR intervals from the other M-L Tang et al paper, Tables 4, 6 & 8
  # (cf methods AS = Score:14 or NB:15, MOVER-W = WCI:3, MOVER-J ≠ JCI:4)
  xs <- rbind(c(8, 3, 1, 2),
              c(22, 2, 0, 1),
              c(43, 0, 1, 0))
  egCI <- allpairci(x = xs, contrast = "RR",
                    methods <- c("AS", "MOVER-W", "MOVER-J", "MOVER-NJ"),
                    alpha=0.05)
  ftable(round(egCI, 4), row.vars = 1, col.vars=c(3, 2))

  # Check qbeta vs qf version of Jeffreys interval, for p = 11/14
  x <- 11
  n <- 14
  alpha <- 0.05
  cc <- 0
  ai <- bi <- 0.5
  CI_qf <- c((2*x + 1) / (2*x + 1 + (2*(n-x) + 1)*qf(1-alpha/2, 2*(n-x)+1, 2*x+1)),
             (2*x + 1) / (2*x + 1 + (2*(n-x) + 1)*qf(alpha/2, 2*(n-x)+1, 2*x+1)))
  CI_beta <- c(qbeta(alpha / 2, x + (ai - cc), n - x + (bi + cc)),
               qbeta(1 - alpha / 2, x + (ai + cc), n - x + (bi - cc)))
  rbind(CI_qf, CI_beta)


  est <- qbeta(0.5, x + (ai), n - x + (bi)) # Obtain phat as the median
  #    CI_upper <- qbeta(1 - alpha / 2, x + (ai + cc), n - x + (bi - cc))
  CI_upper <- (2*x + 1) / (2*x + 1 + (2*(n-x) + 1)*qf(alpha/2, 2*(n-x)+1, 2*x+1))



  #############################################################################
  ### p.12 evaluation of type I error rate (TIER) for test for association
  #############################################################################

  # Load the saved dataset
  load(file=paste0(outpath,"tiers.Rdata"))

  # TIER summaries matching Fagerland 2013
  apply(mytiers[,4:5], 2, mean)
  apply(mytiers[,4:5], 2, max)
  apply(mytiers[,4:5], 2, function(x) mean(x > 0.05))
  apply(mytiers[,4:5], 2, function(x) mean(x < 0.03))




  ### OPTIONAL: run the code below to reproduce the analysis,
  ### or run with different set of parameters

  # PLACEHOLDER: simplify to use test formula from paper instead of scorepairci()

  tiern <- function(ns = fixn,
                    myparams = myparams) {

    tierout <- NULL
    n <- 10
    for (n in ns) {
      cat(paste0("N=", n,"\n"))

    g <- expand.grid(x11 = 0:n, x12 = 0:n, x21 = 0:n)
    # reduce to possible combinations of a,b,c for paired data.
    g <- g[(g$x12 <= n - g$x11) &
             (g$x21 <= n - g$x11 - g$x12), ]
    xs <- data.matrix(cbind(g, x22 = n - rowSums(g)))
    ndis <- rowSums(xs[, 2:3])

    psi <- phi <- NULL
    pbapply::pboptions(style=1)
    i <- 1
    out <- pbapply::pbsapply(1:dim(myparams)[1], function(i)   {
      if (dimnames(myparams)[[2]][2] == "psi") {
        psi <- myparams[i, 2]
      }
      if (dimnames(myparams)[[2]][2] == "phi") {
        phi <- myparams[i, 2]
      }
#      psi <- myparams[i, 2]
      p1 <- p2 <- myparams[i, 1]
      prob <- pdfpair(p1 = p1,
                      p2 = p2,
                      psi = psi,
                      phi = phi,
                      x = xs)
#      sel <-
      xsub <- xs[prob > 1E-8, , drop=F]
      ndissub <- ndis[prob > 1E-8]
      #    dim(xsub)
      if (dim(xsub)[1] > 0) {

        # vectorised 'N-1' test from equation
        n1test <- pchisq(((n-1)/n) * (xsub[, 3] - xsub[, 2])^2 / ndissub, df = 1, lower.tail = FALSE)
        n1test[ndissub == 0] <- 1
        tier <- (n1test < 0.05) %*% prob[prob > 1E-8]

        # vectorised McNemar test from equation
        mactest <- pchisq((xsub[, 3] - xsub[, 2])^2 / ndissub, df = 1, lower.tail = FALSE)
        mactest[ndissub == 0] <- 1
        tier2 <- (mactest < 0.05) %*% prob[prob > 1E-8]
if (FALSE) {
        pvals2 <- sapply(1:dim(xsub)[[1]], function(i)
          pchisq(scorepair(theta = 0,
                           x = xsub[i,],
                           contrast = "RD",
                           cc = FALSE,
                           skew = FALSE,
                           bcf = FALSE)$score^2, df=1, lower.tail=F)
        )
        tier2a <- (pvals2 < 0.05) %*% prob[prob > 1E-8]
}

        # McNemar mid-p test from Fagerland 2013
        px <- 2 * pbinom(pmin(xsub[, 2], xsub[, 3]), ndissub, 0.5, lower.tail = TRUE)
#        px <- 2 * pbinom(apply(xsub[, 2:3], 1, min), ndissub, 0.5, lower.tail = TRUE)
        midp <-  px - dbinom(xsub[, 2], ndissub, 0.5)
        midp[xsub[, 2] == xsub[, 3]] <- (1 - 0.5*dbinom(xsub[, 2], ndissub, 0.5))[xsub[, 2]==xsub[, 3]]
        tier3 <- (midp < 0.05) %*% prob[prob > 1E-8]

        # vectorised cc'd McNemar test from equation
        ccmactest <- pchisq((abs(xsub[, 3] - xsub[, 2]) - 1)^2 / ndissub, df = 1, lower.tail = FALSE)
        ccmactest[ndissub == 0] <- 1
        tier4 <- (ccmactest < 0.05) %*% prob[prob > 1E-8]

        # Exact unconditional test doesnt cope with b=c=0
        # and takes too long anyway
      if (FALSE) {
        pvals.exact <- sapply(1:dim(xsub)[[1]], function(j) {
          contingencytables::McNemar_exact_unconditional_test_paired_2x2(matrix(c(xsub[1, ]), nrow=2))$Pvalue
        })
        tier5 <- (pvals.exact < 0.05) %*% prob[prob > 1E-8]
      }
#        summary(midp[xsub[, 2] > xsub[, 3]])
#        summary(midp[xsub[, 2] < xsub[, 3]])
#        summary(midp)

      } else {
        tier <- 0
        tier2 <- 0
        tier2a <- 0
        tier3 <- 0
        tier4 <- 0
        tier5 <- 0
      }

    c(nminus1 = tier, mcnemar = tier2, midp = tier3, mcnemarcc = tier4)
    })

    tierout <- rbind(tierout, cbind(n = n, myparams, t(out)))
    }

    tierout
  }

  myparams <- expand.grid(p1 = 0.2, psi = 3)
  myparams <- expand.grid(p1 = 0.2, phi = 0.25)

  # Parameter scenarios matching Fagerland 2013
  myparams1 <- expand.grid(p1 = seq(0, 1, 0.01), psi = c(1, 2, 3, 5, 10))
  system.time(tiers1 <- tiern(ns = seq(10, 100, 5), myparams = myparams1))[[3]]/60
  # TIER summaries matching Fagerland 2013
  # Note 'N-1' test performance is almost identical to the exact unconditional method
  apply(tiers1[,4:7], 2, mean)
  apply(tiers1[,4:7], 2, max)
  apply(tiers1[,4:7], 2, function(x) mean(x > 0.05))
  apply(tiers1[,4:7], 2, function(x) mean(x < 0.03))

  # Extended parameter combinations with stronger correlations
  myparams2 <- expand.grid(p1 = seq(0, 1, 0.02), phi = seq(0.25, 0.75, 0.05))
  system.time(tiers2 <- tiern(ns = seq(10, 100, 5), myparams = myparams2))[[3]]/60

  # Combined extended parameter combinations (Note: negative values of phi lead to NAs)
  # Runtime: 46 mins
  myparams3 <- expand.grid(p1 = seq(0, 1, 0.02), phi = c(0, seq(0.05, 0.95, 0.1)))
  nseq <- seq(10, 200, 10)
  nseq2 <- nseq + floor(runif(length(nseq),-4, 6))
  system.time(tiers3 <- tiern(ns = nseq2, myparams = myparams3))[[3]]/60

  # Unexplained issue with one parameter combination needs checking:
 #  n  p1  phi nminus1 mcnemar midp mcnemarcc
 #143 0.5 0.25       0       0    0         0


  save(tiers1, file = paste0(outpath, "newtiers1.Rdata"))
  save(tiers2, file = paste0(outpath, "newtiers2.Rdata"))
  save(tiers3, file = paste0(outpath, "newtiers3.Rdata"))


  # Create a plot of TIERs
  mytiers <- tiers3
  res.factor <- 3
  grid.factor <- 2
  tiff(file = paste0(outpath,"_tiff/","Laud_Fig5new3.tiff"),
       width = 300*grid.factor*res.factor,
       height = 600*res.factor,
       type = "windows"
       #       type="quartz"
  )
  #  par(pty='s')
  par(mfrow = c(2, 2))
  par(cex.main = grid.factor*res.factor*0.8*1, cex.axis=grid.factor*res.factor*0.5*1)
  #  par(mar = res.factor*(c(2,3,3,0.5)+0.1))
  methods <- c("nminus1", "midp", "mcnemar", "mcnemarcc")
  labels <- c("\'N - 1\' AS", "mid-p", "McNemar asymptotic", "McNemar asymptotic (cc)")
  for (i in 1:4) {
    par(mar = grid.factor*res.factor*(c(2,3,3,0.5)+0.1))
    plot(mytiers$p1,
         eval(parse(text=paste0("mytiers$", methods[i]))),
         type = "n",
         ylim = c(0, 0.06),
         xlab = '',
         ylab = '',
         main = labels[i],
         xaxt='n',
         yaxt='n',
         cex.lab = res.factor
    )
    axis(side = 2, las = 2)
    axis(side = 1, las = 1, padj=1)
    mtext(side = 1,
          text = bquote(paste(italic(p)[1]," = ",italic(p)[2])),
          cex = res.factor*1,
          line = 1.5*1.5*res.factor)
    mtext(side = 2,
          text = "Type I error rate",
          cex = res.factor*1,
          line = 1.5*2*res.factor)
    abline(h=0.05, lty=3, lwd=res.factor)

    #    for (ps in c(1, 2, 3, 5, 10)) {
    for (ps in unique(mytiers[,3])) {
      for (n in nseq2[1:4]) {
        #        tiersub <- mytiers[mytiers$psi == ps & mytiers$n == n, ]
        tiersub <- mytiers[mytiers[,3] == ps & mytiers$n == n, ]
        lines(tiersub$p1,
              eval(parse(text=paste0("tiersub$", methods[i]))),
              lty = 2,
              lwd = 0.5*res.factor,
              col = "gray50")
      }
      for (n in nseq2[5:length(nseq2)]) {
        #        tiersub <- mytiers[mytiers$psi == ps & mytiers$n == n, ]
        tiersub <- mytiers[mytiers[,3] == ps & mytiers$n == n, ]
        lines(tiersub$p1,
              eval(parse(text=paste0("tiersub$", methods[i]))),
              lty = 1,
              lwd = 0.25*res.factor)
      }
    }

  }
  dev.off()




  #############################################################################
  ### SUPPLEMENTARY FIGURES:
  ### CP, MACP, location index RNCP and width for selected methods for RD,
  ### with N = 40, \alpha=0.05 and \phi=0.25
  #############################################################################
  # RD
  RDpairteam <- c("SCAS-bc", "AS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP", "AM", "Wald") #Paired RD
  RDcpairteam <- c("SCAS-c125", "SCAS-c25", "SCAS-c5", "MOVER-c125", "MOVER-c25", "MOVER-c5", "Wald-cc") 	#Paired RD, cc
  teamlist <- list(RDpairteam, RDcpairteam)
  teamlabels <- c("RDpair", "RDcpair")
#  load(file = paste0(outpath, "cparrays.RD.", 40, ".", 200, ".Rdata"))
  for (n in c(20, 40, 65)) {
  load(file = paste0(outpath, "cparrays.RD.", n, ".", 200, ".Rdata"))
  for (j in c(0.1, 0.25, 0.5, 0.75)) {
#    for (j in c(-0.25, -0.1, 0.9, 0.99)) {
    for (i in c(0.05, 0.1, 0.01)) {
        #  for (i in c(0.05)) {
        for (k in 1:2) {
          if (!((k ==2) & (n %in% c(20, 65)))) {
          plotpanel(
            plotdata = arrays, alpha = i, par3 = j,
            limits = c(0, 1), sel = teamlist[[k]], plotlab = teamlabels[k],
            fmt = "png", res.factor = 4, CIlen = TRUE
          )
          }
        }
      }
    }
  }

  # selection of RR methods for supplementary plots
  RRpairteam <- c("SCAS-bc", "AS-bc", "AS", "MOVER-NJ", "MOVER-NW", "MOVER-W", "BP", "BP-J", "Wald")  	#Paired RR
  RRcpairteam <- c("SCAS-c125", "SCAS-c25", "SCAS-c5", "MOVER-c125", "MOVER-c25", "MOVER-c5", "BP-cc") 	#Paired RR, cc
  teamlist <- list(RRpairteam, RRcpairteam)
  teamlabels <- c("RRpair", "RRcpair")
#  load(file=paste0(outpath, "cparrays.RR.", 40, ".",200,".Rdata"))
  for (n in c(20, 40, 65)) {
    load(file=paste0(outpath, "cparrays.RR.", n, ".",200,".Rdata"))
    for (j in c(0.1, 0.25, 0.5, 0.75)) {
      for (i in c(0.05, 0.1, 0.01)) {
        for (k in 1:2) {
          if (!((k ==2) & (n %in% c(20, 65)))) {
            plotpanel(plotdata=arrays, alpha=i, par3=j,
                    limits=c(0,1), sel=teamlist[[k]], plotlab=teamlabels[k],
                    fmt="png", res.factor = 4, CIlen = TRUE)
          }
        }
      }
    }
  }

  # OR
  ORpairteam <- c("SCASp", "SCASpu", "mid-p", "Jeffreys", "Wilson", "Wald", "Laplace")
  ORcpairteam <- c("SCASp-c125", "SCASp-c25", "SCASp-c5", "midp-c25", "Jeffreys-c125", "Jeffreys-c25", "C-P", "Blaker")
  teamlist <- list(ORpairteam, ORcpairteam)
  teamlabels <- c("ORpair", "ORcpair")
#  load(file=paste0(outpath, "cparrays.OR.", 40, ".",200,".Rdata"))
  for (n in c(20, 40, 65)) {
    load(file=paste0(outpath, "cparrays.OR.", n, ".",200,".Rdata"))
    for (j in c(0.1, 0.25, 0.5, 0.75)) {
      for (i in c(0.05, 0.1, 0.01)) {
        for (k in 1:2) {
          if (!((k ==2) & (n %in% c(20, 65)))) {
            plotpanel(plotdata=arrays, alpha=i, par3=j,
                  limits=c(0,1), sel=teamlist[[k]], plotlab=teamlabels[k],
                  fmt="png", res.factor = 4, CIlen = TRUE)
          }
        }
      }
    }
  }


  # Sample code to retrieve a previously run array:
  # load(file=paste(outpath, "cis.OR.", 105, ".Rdata",sep=""))
  # mycis <- arrays
  # load(file=paste(outpath, "cparrays.RR.", 40, ".",200,".Rdata",sep=""))
  # arr <- myarrays



  if(FALSE) {
    # Sample code in case needed to combine plots for publication
    # install.packages("png")
    library(png)
    fig1a <- readPNG(paste0(outpath,"_png/summaryRD95_150,50os.png"), TRUE)
    fig1b <- readPNG(paste0(outpath,"_png/summaryID95_150,50os.png"), TRUE)
    res.factor <- 6
    tiff(file=paste0(outpath,"_tiff/Figure1.tiff"),width=600*res.factor,height=2*330*res.factor,type="quartz")
    # png(file=paste0(outpath,"_png/Figure1.png"),width=600*res.factor,height=2*330*res.factor,type="quartz")
    par(mar=c(0,0,0,0),oma=c(0,0,0,0))
    plot(0:1,0:1, type='n',axes=F, xaxs="i", yaxs="i")
    # plot(0:1,0:1, axes=F, xaxs="i", yaxs="i")
    rasterImage(fig1a, 0, 0.5, 1, 1)
    rasterImage(fig1b, 0, 0, 1, 0.5)
    text(0.02,0.98,"(a)",cex=1.5*res.factor)
    text(0.02,0.48,"(b)",cex=1.5*res.factor)
    dev.off()
  }


}
