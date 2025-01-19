
set.seed(2012) #ensure we use the same jitters for each run


RRpairteam <- c("Tang", "Tang-bc", "Tang-sc", "Tang-scbc") 	#Paired RR
RR2pairteam <- c("MOVER-w", "MOVER-nw", "MOVER-nj","MOVER-ns") 	#Paired RR
RRccpairteam <- c("Tang-cc5", "Tang-cc125", "Tang-sccc", "MOVER-ccns") 	#Paired RR, cc

RDpairteam <- c("Tango", "Tango-bc", "Tango-sc", "Tango-scbc") 	#Paired RD
RD2pairteam <- c("MOVER-w", "MOVER-nw", "MOVER-nj","MOVER-ns") 	#Paired RD
RDccpairteam <- c("Tango-cc5", "Tango-cc125", "Tango-sccc", "MOVER-ccns") 	#Paired RD, cc

#Remotes:
#  petelaud/ratesci-dev

#source("/Users/ssu/Documents/Main/GitHub/ratesci/R/pairbinci.R")
#source("/Users/ssu/Documents/Main/GitHub/ratesci/R/moverci.R")
#source("/Users/ssu/Documents/Main/GitHub/ratesci/R/scoreci.R")
#source("/Users/ssu/Documents/Main/GitHub/ratesci/R/rateci.R")

root <- "/Users/ssu/Documents/"
outpath <- paste(root, "Main/Courses_papers/skewscore/paired/", sep = "")

if (FALSE) {

  #### RD
  # n=20
  # 1.5min for cis per alpha incl skew
  # +6min per combo with n=200 with smoothing
  # n=40
  # 11min per alpha (33)
  # +39min per combo for n=200 with smoothing (8h for 3*4)

#system.time(mycis <- cifun(n=10, contrast="RR", n.grid=100,
                                 #alph = c(0.01, 0.05, 0.1), psis = c(1, 2, 10, 100)))[[3]]/60
#                                 alph = c(0.01, 0.05, 0.1), psis = c(1, 2, 10, 100)))[[3]]/60
system.time(mycis <- cifun(n=40, contrast="RR", alph = c(0.05, 0.1)))[[3]]/60
Sys.time()
system.time(
  cparrays_RR40 <- cpfun(ciarrays = mycis,
                         n.grid=200, phis=c(0, 0.25, 0.5), smooth=FALSE)
  )[[3]]/60

teamlist <- list(RDpairteam, RD2pairteam, RDccpairteam)
teamlabels <- c("RDpair", "RD2pair", "RDccpair")

teamlist <- list(RRpairteam, RR2pairteam, RRccpairteam)
teamlabels <- c("RRpair", "RR2pair", "RRccpair")

# load(file=paste0(outpath, "cparrays.RR.", 40, ".",200,".Rdata"))
# cparrays_RR40 <- arrays
# dimnames(cparrays_RR40)[[1]]
for (j in c(0, 0.25, 0.5)) {
  for (i in c(0.05, 0.1)) {
    for (k in 1:2) {
      plotpanel(plotdata=cparrays_RR40, alpha=i, par3=j,
                limits=c(0,1), sel=teamlist[[k]], oneside=F, plotlab=teamlabels[k],
                res.factor=6, fmt="png", linesx=F, colour=T, sided="R",
                smoothed=FALSE)
    }
  }
}


x <- c(1,1,7,12)
round(pairbinci(x=x, contrast = "RR", method_RR = "Score_closed")$estimates[,c(1,3)], 3)
round(pairbinci(x=x, contrast = "RR", method_RR = "Score", skew = T)$estimates[,c(1,3)], 3)
round(pairbinci(x=x, contrast = "RR", method_RR = "MOVER", moverbase = "SCAS")$estimates[,c(1,3)], 3)
round(pairbinci(x=x, contrast = "RR", method_RR = "MOVER", moverbase = "jeff")$estimates[,c(1,3)], 3)




}
