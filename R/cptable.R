
if (FALSE) {

  # Combine output arrays for summarising across different Ns
  mynums <- c(20, 40, 65) #, 65, 105)
  mymethods <- c("SCAS-bc", "SCAS", "AS", "MOVER-NJ", "MOVER-W", "BP")
#  submethods <- c("SCAS-bc", "AS", "MOVER-NJ", "MOVER-W", "BP")
  nmeth <- length(mymethods)
  load(file=paste0(outpath, "cparrays.RR.", 20, ".",200,".Rdata"))
  #  dimnames(arrays$summaries)[[2]][16] <- "BP"
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
#  dimnames(arrays$summaries)

num <- 40
  submethods <- mymethods
  for (num in c(20, 40)) {
    load(file=paste0(outpath, "cparrays.RD.", num, ".",200,".Rdata"))
    bigarray[,submethods,,,paste(num), "RD"] <- arrays$summaries[,submethods,,,paste(num),]
    load(file=paste0(outpath, "cparrays.RR.", num, ".",200,".Rdata"))
    dimnames(arrays$summaries)[[2]][16] <- "BP"
    bigarray[,submethods,,,paste(num), "RR"] <- arrays$summaries[,submethods,,,paste(num),]
  }
#  submethods <- mymethods[-2]
  num <- "65"
  for (num in "65") {
    load(file=paste0(outpath, "cparrays.RD.", num, ".",200,".Rdata"))
    bigarray[,submethods,,,paste(num), "RD"] <- arrays$summaries[,submethods,,,paste(num),]
    load(file=paste0(outpath, "cparrays.RR.", num, ".",200,".Rdata"))
    bigarray[,submethods,,,paste(num), "RR"] <- arrays$summaries[,submethods,,,paste(num),]
  }

dimnames(arrays$summaries)

#library(xtable)
### Tables for publication
#outpath=paste(root,"Main/Courses_papers/skewscore/plots/",sep="")
#outorder<-c("SC","MN","MOVER-J","Wald")
#selectmethods <- c("SCAS", "SCAS-bc", "AS", "MOVER-NJ", "MOVER-W", "BP-W")
#selectmethods <- c("SCAS", "SCAS-bc", "AS", "MOVER-NJ", "MOVER-W", "BP")
selectalpha <- c("95", "90", "99")
#numslist<-c("20","40")
#numslist<-c("20")

dimnames(bigarray)

#load(file=paste(outpath,"masterarraysx.",m,".",n,".",n.grid,".Rdata",sep=""))
#load(file=paste0(outpath, "cparrays.RD.", 40, ".",200,".Rdata"))
#load(file=paste0(outpath, "cparrays.RR.", 40, ".",200,".Rdata"))

#dimnames(arrays$summaries)
mysummaries <-
  bigarray[, selectmethods, selectalpha, c("meanCP", "pctCons", "pctnear", "pctgoodloc", "meanlocindex"),,]
#dim(mysummary)

#mysummaries[,,,,"65","RD"]

mysummaries[, , , c("pctCons", "pctnear", "pctgoodloc"),,] <- round(as.numeric(mysummaries[, , ,c("pctCons", "pctnear", "pctgoodloc"),,]), 0)

anticons <- (as.numeric(mysummaries[,,,"pctCons",,]) < 50)
mysummaries[,,,"pctnear",,] <- paste0(ifelse(anticons,"["," "),
                                 mysummaries[,,,"pctnear",,],
                                 ifelse(anticons,"]"," "))


# Create flat table from multidimensional array
mytable <- ftable((mysummaries[, , , c("pctnear", "pctgoodloc")]), col.vars = c(1,4), row.vars = c(3,2))
mytable1 <-
  ftable((mysummaries[, , , c("pctnear"),,]), col.vars = c(3,1), row.vars = c(5,4,2))
write.ftable(mytable1, sep=',', quote=TRUE, file = paste0(outpath, "tablecp2065.csv"))

mytable2 <-
  ftable((mysummaries[, , , c("pctgoodloc"),,]), col.vars = c(3,1), row.vars = c(5,4,2))
write.ftable(mytable2, sep=',', quote=TRUE, file = paste0(outpath, "tableloc2065.csv"))



# MACP table

mysummaries <-
  (arrays$summaries)[, selectmethods, ,c("pctCons", "pctnear", "pctAvenear", "pctgoodloc", "meanlocindex"),,]
#dim(mysummary)

mysummaries[, , ,c("pctCons", "pctnear", "pctAvenear", "pctgoodloc")] <- round(as.numeric(mysummaries[, , ,c("pctCons", "pctnear", "pctAvenear", "pctgoodloc")]), 0)

anticons <- (as.numeric(mysummaries[,,,"pctCons"]) < 50)
mysummaries[,,,"pctAvenear"] <- paste0(ifelse(anticons,"["," "),
                                       mysummaries[,,,"pctAvenear"],
                                       ifelse(anticons,"]"," "))
mysummaries[,,,"pctnear"] <- paste0(ifelse(anticons,"["," "),
                                       mysummaries[,,,"pctnear"],
                                       ifelse(anticons,"]"," "))


# Create flat table from multidimensional array
for (summ in c("pctnear", "pctAvenear", "pctgoodloc", "meanlocindex")) {
  mytable <- ftable((mysummaries[, , , summ]), col.vars = c(1), row.vars = c(3,2))
  write.ftable(mytable, sep=',', quote=TRUE, file = paste0(outpath, "tableRD40",summ,".csv"))
  mytable
}



load(file=paste(outpath,"bignsummary.Rdata"))
mytable2 <- ftable(bignsummary[,"dncp",,,,], col.vars = c(3,2), row.vars = c(4,1))
write.ftable(mytable2, sep=',', quote=TRUE, file = paste0(outpath, "bigntable205.csv"))







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

}
