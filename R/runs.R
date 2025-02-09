
set.seed(2012) #ensure we use the same jitters for each run


RRpairteam <- c("SCAS-bc", "SCAS", "AS", "MOVER-NJ", "MOVER-W", "BP") 	#Paired RR
#RRpairteam <- c("SCAS", "SCAS-bc", "AS", "AS-bc") 	#Paired RR
#RRpairteam <- c("SCAS", "AS", "MOVER-NJ", "MOVER-W") 	#Paired RR
#RR2pairteam <- c("MOVER-NJ", "MOVER-W", "BP-J", "BP-W") 	#Paired RR
#RR2pairteam <- c("MOVER-w", "MOVER-nw", "MOVER-nj","MOVER-ns") 	#Paired RR
RRccpairteam <- c("SCAS-cc125", "SCAS-cc5", "MOVER-NJcc125", "MOVER-NJcc5") 	#Paired RR, cc

#RDpairteam <- c("SCAS-bc", "AS", "MOVER-NJ", "MOVER-W", "BP") 	#Paired RD
RDpairteam <- c("SCAS-bc", "SCAS", "AS", "MOVER-NJ", "MOVER-W", "BP") 	#Paired RD
#RDpairteam <- c("SCAS", "SCAS-bc", "AS", "AS-bc" ) 	#Paired RD
#RDpairteam <- c("SCAS", "Tango", "MOVER-nj", "MOVER-w") 	#Paired RD
#RD2pairteam <- c("MOVER-NJ","MOVER-W", "MOVER-NW", "BP") 	#Paired RD
#RDccpairteam <- c("Tango-cc5", "Tango-cc125", "SCAS-cc", "MOVER-NJcc") 	#Paired RD, cc
RDccpairteam <- c("SCAS-cc125", "SCAS-cc5", "MOVER-NJcc125", "MOVER-NJcc5") 	#Paired RD, cc

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

load(file=paste0(outpath, "cis.RD.", 41, ".Rdata"))
mycis <- ciarrays
#system.time(mycis <- cifun(n=10, contrast="RR", n.grid=100,
                                 #alph = c(0.01, 0.05, 0.1), psis = c(1, 2, 10, 100)))[[3]]/60
#                                 alph = c(0.01, 0.05, 0.1), psis = c(1, 2, 10, 100)))[[3]]/60
system.time(mycis <- cifun(n=40, contrast="RD", alph = c(0.1, 0.05, 0.01)))[[3]]/60
#dimnames(mycis[[2]])
Sys.time()
system.time(
  arrays <- cpfun(ciarrays = mycis,
                      n.grid=200, phis=c(0.1, 0.25, 0.5, 0.75), smooth=TRUE)
          )[[3]]/60

teamlist <- list(RDpairteam,  RD2pairteam)
teamlabels <- c("RDpair", "RD2pair")

teamlist <- list(RRpairteam, RR2pairteam)
teamlabels <- c("RRpair", "RR2pair")

teamlist <- list(RDpairteam)
teamlabels <- c("RDpair")

teamlist <- list(RRpairteam, RRccpairteam)
teamlabels <- c("RRpair", "RRccpair")

teamlist <- list(RDpairteam, RDccpairteam)
teamlabels <- c("RDpair", "RDccpair")

# N=40 plugged in: RD 3 mins per alpha with 5 methods, xx mins per phi with 200grid +smooth
# Larger N
# N=65 plugged in: RD 4 mins per alpha with (wrong) 5 methods, 41 mins per phi with 200grid +smooth
# N=65 on battery: RD 12 mins per alpha with 6 methods, 45 min per phi
# N=65 plugged in: RR 11 mins per alpha with 5 methods, 45 mins per phi with 200grid +smooth
# N=105 46 mins per alpha with 5 methods

system.time(mycis <- cifun(n=10, contrast="RD", alph = alphas, methods = RDpairteam))[[3]]/60
mycis$cis[1:10,, "SCAS","95", "10", "RD"]
dimnames(mycis$cis)


RDpairteam <- RRpairteam <- c("SCAS-bc", "SCAS", "AS", "MOVER-NJ", "MOVER-W", "BP") 	#Paired RR
alphas <- c(0.1, 0.05, 0.01)
phis <- c(0.1, 0.25, 0.5, 0.75)
if (FALSE) {
  system.time(mycis <- cifun(n=10, contrast="RD", alph = 0.05, methods = RDpairteam))[[3]]/60
  system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=0.25))[[3]]/60
  system.time(mycis <- cifun(n=10, contrast="RR", alph = 0.05, methods = RRpairteam))[[3]]/60
  system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=0.25))[[3]]/60
  system.time(mycis <- cifun(n=40, contrast="RD", alph = 0.05, methods = RDpairteam))[[3]]/60
  system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=0.25))[[3]]/60
  system.time(mycis <- cifun(n=40, contrast="RR", alph = 0.05, methods = RRpairteam))[[3]]/60
  system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=0.25))[[3]]/60
}
#load(file=paste0(outpath, "cis.RD.", 20,".Rdata"))
#mycis <- ciarrays
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
system.time(mycis <- cifun(n=105, contrast="RD", alph = 0.05, methods = RDpairteam))[[3]]/60
# load(file=paste0(outpath, "cis.RD.", 105,".Rdata"))
#mycis <- ciarrays
# 9 hours
Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60
system.time(mycis <- cifun(n=105, contrast="RR", alph = 0.05, methods = RRpairteam))[[3]]/60
Sys.time(); system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=phis))[[3]]/60

#Sys.time()
#system.time(arrays <- cpfun(ciarrays = mycis, n.grid=200, phis=c(0.25), smooth=TRUE))[[3]]/60
teamlist <- list(RRpairteam); teamlabels <- c("RRpair")
teamlist <- list(RDpairteam); teamlabels <- c("RDpair")
plotpanel(plotdata=arrays, alpha=0.05, par3=0.25,
          sel=teamlist[[1]], plotlab=teamlabels[1],
          fmt="png")
dev.off()
#mycis$cis[1,,,,,]
#allpairci(x=c(1,1,7,12), methods=RDpairteam)

# load(file=paste0(outpath, "cparrays.RR.", 40, ".",200,".Rdata"))
# load(file=paste0(outpath, "cis.RD.", 105,".Rdata"))
dimnames(arrays$mastercp)
mycis$cis[,,"SCAS","95",,]
arrays$mastercp[100:110,100:110,"0.25","SCAS","95","cp",,]
dimnames(mycis$cis)
# cparrays_RD30 <- arrays
# names(cparrays_RD30)
# dimnames(cparrays_RD30)[3]
# dim(cparrays_RD30)[[3]]
#for (j in c(0.1, 0.25, 0.5)) {
#  for (i in c(0.01, 0.05, 0.1)) {

#for (j in c(0.1, 0.25, 0.5, 0.75)) {
for (j in c(0.1, 0.25, 0.5, 0.75)) {
  for (i in c(0.05, 0.1, 0.01)) {
#    for (j in c(0.25)) {
#      for (i in c(0.05)) {
        for (k in 1) {
      plotpanel(plotdata=arrays, alpha=i, par3=j,
                limits=c(0,1), sel=teamlist[[k]], oneside=F, plotlab=teamlabels[k],
                res.factor=6, fmt="png", linesx=F, colour=T, sided="R",
                smoothed=FALSE)
    }
  }
}

install.packages('latex2exp')
dev.off()

load(file=paste0(outpath, "cis.RD.", 40,".Rdata"))
mycis <- ciarrays
load(file=paste0(outpath, "cparrays.RD.", 40, ".",200,".Rdata"))
#cparrays_RD30 <- arrays
load(file=paste0(outpath, "cparrays.RR.", 40, ".",200,".Rdata"))
#cparrays_RR30 <- arrays$summaries[,-12,,,,]

dim(cparrays_RD30$summaries)
dim(cparrays_RR30$summaries)
dimnames(cparrays_RD30$summaries)[[2]]
dimnames(cparrays_RR30$summaries)[[2]]

mydims <- dim(cparrays_RD30$summaries)
mydims[6] <- 2
mydimnames <- dimnames(cparrays_RD30$summaries)
mydimnames[[6]] <- c("RD", "RR")
bigarray <- array(NA, dim = mydims)
dimnames(bigarray) <- mydimnames
bigarray[,,,,,1] <- cparrays_RD30$summaries
bigarray[,,,,,2] <- cparrays_RR30$summaries[,-12,,,,]

selectmethods <- c("SCAS", "AS", "MOVER-NJ", "MOVER-W")
mysummary <- bigarray[, selectmethods, "95", c("pctnear", "pctgoodloc"),,]

library(data.table)
as.data.table(mysummary)

#𝜙 𝜙 ϕ

install.packages("Unicode")
#φ
special_char <- "ϕ"
Unicode::as.u_char(utf8ToInt(special_char))


par(pty = 's')
for (contrast in c("RD", "RR")) {
for (phi in c("0.1", "0.25", "0.5")) {
  plot(as.numeric(mysummary[phi,,"pctnear", contrast]),
       as.numeric(mysummary[phi,,"pctgoodloc", contrast]),
       xlim = c(0, 100), ylim=c(0, 100),
       xlab = "% proximate",
       ylab = "% central",
       main = paste0(contrast, ", \u03D5=", phi)
       )
  text(as.numeric(mysummary[phi,,"pctnear", contrast]),
       as.numeric(mysummary[phi,,"pctgoodloc", contrast]),
       labels = dimnames(mysummary)[[2]],
       pos = 4
       )
}}

# Large sample size examples
onecpfun(0.3, 0.1, n=30, contrast = "RD", alph=0.05, phis=0.25, methods=RRpairteam)
system.time(cp205RD25 <- onecpfun(0.4, 0.1, n=205, contrast = "RD", alph=0.05, phis=0.25, methods=RDpairteam))[[3]]/60
system.time(cp205RD75 <- onecpfun(0.3, 0.2, n=205, contrast = "RD", alph=0.05, phis=0.75, methods=RDpairteam))[[3]]/60
system.time(cp205RD2599 <- onecpfun(0.4, 0.1, n=205, contrast = "RD", alph=0.01, phis=0.25, methods=RDpairteam))[[3]]/60
system.time(cp205RD7599 <- onecpfun(0.3, 0.2, n=205, contrast = "RD", alph=0.01, phis=0.75, methods=RDpairteam))[[3]]/60
system.time(cp205RD10 <- onecpfun(0.4, 0.1, n=205, contrast = "RD", alph=0.05, phis=0.10, methods=RDpairteam))[[3]]/60
system.time(cp205RD50 <- onecpfun(0.4, 0.2, n=205, contrast = "RD", alph=0.05, phis=0.50, methods=RDpairteam))[[3]]/60
system.time(cp205RD1099 <- onecpfun(0.4, 0.1, n=205, contrast = "RD", alph=0.01, phis=0.10, methods=RDpairteam))[[3]]/60
system.time(cp205RD5099 <- onecpfun(0.4, 0.2, n=205, contrast = "RD", alph=0.01, phis=0.50, methods=RDpairteam))[[3]]/60

round(cp205RD25[,c(1, 4, 5),], 3)
round(cp205RD75[,c(1, 4, 5),], 3)
round(cp205RD2599[,c(1, 4, 5),], 3)
round(cp205RD7599[,c(1, 4, 5),], 3)
round(cp205RD10[,c(1, 4, 5),], 3)
round(cp205RD50[,c(1, 4, 5),], 3)
round(cp205RD1099[,c(1, 4, 5),], 3)
round(cp205RD5099[,c(1, 4, 5),], 3)

system.time(cp205RR25 <- onecpfun(0.4, 0.1, n=205, contrast = "RR", alph=0.05, phis=0.25, methods=RRpairteam))[[3]]/60
system.time(cp205RR75 <- onecpfun(0.3, 0.2, n=205, contrast = "RR", alph=0.05, phis=0.75, methods=RRpairteam))[[3]]/60
system.time(cp205RR2599 <- onecpfun(0.4, 0.1, n=205, contrast = "RR", alph=0.01, phis=0.25, methods=RRpairteam))[[3]]/60
system.time(cp205RR7599 <- onecpfun(0.3, 0.2, n=205, contrast = "RR", alph=0.01, phis=0.75, methods=RRpairteam))[[3]]/60
system.time(cp205RR10 <- onecpfun(0.4, 0.1, n=205, contrast = "RR", alph=0.05, phis=0.10, methods=RDpairteam))[[3]]/60
system.time(cp205RR50 <- onecpfun(0.4, 0.2, n=205, contrast = "RR", alph=0.05, phis=0.50, methods=RDpairteam))[[3]]/60
system.time(cp205RR1099 <- onecpfun(0.4, 0.1, n=205, contrast = "RR", alph=0.01, phis=0.10, methods=RDpairteam))[[3]]/60
system.time(cp205RR5099 <- onecpfun(0.4, 0.2, n=205, contrast = "RR", alph=0.01, phis=0.50, methods=RDpairteam))[[3]]/60

round(cp205RR25[,c(1, 4, 5),], 3)
round(cp205RR75[,c(1, 4, 5),], 3)
round(cp205RR2599[,c(1, 4, 5),], 3)
round(cp205RR7599[,c(1, 4, 5),], 3)
round(cp205RR10[,c(1, 4, 5),], 3)
round(cp205RR50[,c(1, 4, 5),], 3)
round(cp205RR1099[,c(1, 4, 5),], 3)
round(cp205RR5099[,c(1, 4, 5),], 3)

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

bignsummary[,"dncp",,"0.05",,]
save(bignsummary, file = paste(outpath, "bignsummary.Rdata"))



round(cbind(
  cp205RR10[,4,],
  cp205RR25[,4,],
  cp205RR50[,4,],
  cp205RR75[,4,],
  cp205RR1099[,4,],
  cp205RR2599[,4,],
  cp205RR5099[,4,],
  cp205RR7599[,4,]
), 4)



x <- c(1,1,7,12)
round(pairbinci(x=x, contrast = "RR", method_RR = "Score_closed")$estimates[,c(1,3)], 3)
round(pairbinci(x=x, contrast = "RR", method_RR = "Score", skew = T)$estimates[,c(1,3)], 3)
round(pairbinci(x=x, contrast = "RR", method_RR = "MOVER", moverbase = "SCAS")$estimates[,c(1,3)], 3)
round(pairbinci(x=x, contrast = "RR", method_RR = "MOVER", moverbase = "jeff")$estimates[,c(1,3)], 3)


# Tang version of Jeffreys
Y <- 9
n <- 14
(2*Y + 1) / (2*Y + 1 + (2*(n-Y)+1) * qf(0.025, 2*(n-Y)+1, 2*Y+1, lower.tail=T ) )
(2*Y + 1) / (2*Y + 1 + (2*(n-Y)+1) * qf(0.975, 2*(n-Y)+1, 2*Y+1, lower.tail=T ) )
jeffreysci(9,14, adj=T)
x <- c(43, 0, 1, 0)
x <- c(8,3,1,2)
x <- c(4,9,3,16)
x <- c(16,9,3,4)
pairbinci(x=x, contrast="RD", method_RD = "Score", skew=F)$estimates[,c(1,3)]
#pairbinci(x=x, contrast="RD", method_RD = "MOVER", moverbase = "wilson")$estimates[,c(1,3)]
pairbinci(x=x, contrast="RD", method_RD = "MOVER_newc", moverbase = "wilson")$estimates[,c(1,3)]
pairbinci(x=x, contrast="RD", method_RD = "MOVER", moverbase = "jeff")
pairbinci(x=x, contrast="RD", method_RD = "MOVER", moverbase = "jeff", cc=0.125)
pairbinci(x=x, contrast="RD", method_RD = "MOVER_newc", moverbase = "jeff")$estimates[,c(1,3)]
pairbinci(x=x, contrast="RD", method_RD = "MOVER_newc", moverbase = "jeff", cc=0.125)$estimates[,c(1,3)]
pairbinci(x=x, contrast="RD", method_RD = "MOVER_newc", moverbase = "jeff", level=0.95)$estimates[,c(1,3)]
moverpair(x=x, contrast="RD", method = "wilson")
moverpair(x=x, contrast="RD", method = "wilson", corc=T)
moverpair(x=x, contrast="RD", method = "jeff")
moverpair(x=x, contrast="RD", method = "jeff", corc=T, crudemid=F)
moverpair(x=x, contrast="RD", method = "jeff", corc=T, crudemid=T)


x <- c(1, 1, 7, 12)


pairbinci(x=x, contrast="RR", method_RR = "Score", skew=F, bcf=F)
pairbinci(x=x, contrast="RR", method_RR = "MOVER", moverbase = "wilson")
pairbinci(x=x, contrast="RR", method_RR = "MOVER", moverbase = "jeff")


x[1] * x[4] - x[2] * x[3]

?FDist


}
