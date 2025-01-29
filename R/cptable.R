
if (FALSE) {


#library(xtable)
### Tables for publication
#outpath=paste(root,"Main/Courses_papers/skewscore/plots/",sep="")
outorder<-c("SC","MN","MOVER-J","Wald")
selectmethods <- c("SCAS", "SCAS-bc", "AS", "MOVER-NJ", "MOVER-W", "BP-W")
selectmethods <- c("SCAS", "SCAS-bc", "AS", "MOVER-NJ", "MOVER-W", "BP")
numslist<-c("20","40")
numslist<-c("20")

#load(file=paste(outpath,"masterarraysx.",m,".",n,".",n.grid,".Rdata",sep=""))
load(file=paste0(outpath, "cparrays.RD.", 40, ".",200,".Rdata"))
load(file=paste0(outpath, "cparrays.RR.", 40, ".",200,".Rdata"))

#dimnames(arrays$summaries)
mysummaries <-
  (arrays$summaries)[, selectmethods, ,c("meanCP", "pctCons", "pctnear", "pctgoodloc", "meanlocindex"),,]
#dim(mysummary)

mysummaries[, , ,c("pctCons", "pctnear", "pctgoodloc")] <- round(as.numeric(mysummaries[, , ,c("pctCons", "pctnear", "pctgoodloc")]), 0)

anticons <- (as.numeric(mysummaries[,,,"pctCons"]) < 50)
mysummaries[,,,"pctnear"] <- paste0(ifelse(anticons,"["," "),
                                 mysummaries[,,,"pctnear"],
                                 ifelse(anticons,"]"," "))


# Create flat table from multidimensional array
mytable <- ftable((mysummaries[, , ,c("pctnear", "pctgoodloc")]), col.vars = c(1,4), row.vars = c(3,2))
write.ftable(mytable, sep=',', quote=TRUE, file = paste0(outpath, "tableRD40.csv"))



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
