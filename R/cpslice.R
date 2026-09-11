# cpslice function to reproduce 2-D plots e.g. in Fagerland et al.

if (FALSE) {

  root <- "/Users/ssu/Documents/"
  outpath <- paste(root, "Main/Courses_papers/skewscore/paired/", sep = "")
  outpath <- "D:/Pete/Documents/GitHub/cpplot/data/"
  outpath <- "D:/Pete/Documents/Research/paired/" # Remove for final upload


res.factor <- 3
tiff(file = paste0(outpath,"_tiff/unluckyN.tiff"),
width = (600) * res.factor,
height = 200 * res.factor,
#    type="quartz"
type="windows"
)
par(cex.main = res.factor*0.8*1, cex.axis=res.factor*0.8*1)
#par(mar = res.factor*(c(2,3,1,0.5)+0.1))


par(mfrow = c(1, 3)) #, cex = res.factor)

for (myN in c(39, 40, 41)) {

#myN <- 41
psi <- 3
load(file=paste0(outpath, "cis.RR.", myN, ".Rdata"))
#system.time(ciarrays <- cifun(n=myN, contrast="RD", alph = 0.05))[[3]]/60

#p0 <- as.numeric(dimnames(arrays$mastercp)[[1]])
del <- 0.2 # * myN/40
#del <- 0.1
#del <- 0
#length(p0)
#p0 <- seq(0, 1-del, length.out=51)
#p2 <- p0

#p2 <- p0[p0 + del >= 0 & p0 + del <= 1]

# 2-D coverage plot, e.g. illustrating how unrepresentative Fig 8.11 is
p2 <- seq(0, 1-del, length.out=101)
p1 <- p2 + del
  cp1 <- onecpfun(
    n = myN,
    p1 = p1,
    p2 = p2,
#    ciarrays = mycis,
    ciarrays = ciarrays,
    alph = 0.05,
    psis = 3
#    phis = 0.25
  )

  if (FALSE) {
  # One point at a time is inefficient
  cp1 <- onecpfun(
    n = 40,
    contrast = "RD",
    p1 = p1[10],
    p2 = p2[10],
    #    ciarrays = arrays,
    alph = 0.05,
    psis = 3
#        phis = 0.25
  )
  }
  plot(p2,
     cp1[,"AS","cp"],
     type = "l",
     lwd = 2,
     ylim = c(0.90, 1),
     ylab = "Coverage Probability",
     xlab = "p2",
     main = paste0("N = ", myN, ", θ = ", del, "±0.005, ψ = ", psi)
  )
abline(h=0.95)
rect(
  xleft = par("usr")[1], xright = par("usr")[2], ybottom = 0.945, ytop = 0.955,
  border = NA, col = adjustcolor("gray", alpha = 0.3)
)
#abline(h=0.945,lty=2)

dels <- del + seq(-0.005,0.005,0.001)[-6]
for(i in 1:length(dels)){
  p1 <- p2 + dels[i]
  cp1 <- onecpfun(
    p1 = p1,
    p2 = p2,
    ciarrays = ciarrays,
    alph = 0.05,
#    phis = 0.25
    psis = psi
  )
  lines(p2,
       cp1[,"AS","cp"],
       lty = 2
  )
}

}
dev.off()

# 2-D Type I error plot
load(file=paste0(outpath, "cparrays.RD.", 65, ".",200,".Rdata"))
p2 <- p1 <- seq(0,1-del,length.out=101)
cp1 <- onecpfun(
  p1 = p1,
  p2 = p2,
  ciarrays = arrays,
  alph = 0.05,
  phis = 0.25
)


#p <- 0.1
# Attempt to loop onecpfun without reading in an array - slow?
#sapply(seq(0.1,0.3,0.1), function(p) onecpfun(contrast="RD", n=10, p1=p, p2=p, alph=0.05, phis=0.25, methods=c("SCAS", "AS")))

par(pty='s')
plot(p2,
     1 - cp1[,"SCAS-bc","cp"],
     type = "l",
     ylim = c(0, 0.06),
     ylab = "Type I error rate",
     xlab = "p2 = p1",
     main = paste0("N = 65, phi = 0.25"),
     yaxt='n'
)
axis(side = 2, las = 2)
abline(h=0.05, lty=3)
lines(p2, 1 - cp1[,"SCAS-bc","cp"], lty=1, lwd = 2)
lines(p2, 1 - cp1[,"AS","cp"], lty=2, lwd = 2)
lines(p2, 1 - cp1[,"TDAS","cp"], lty=3, lwd = 2)
lines(p2, 1 - cp1[,"SCASstrat","cp"], lty=4, lwd = 2)

# 2-D interval width plot, RD
myN <- 40
load(file=paste0(outpath, "cparrays.RD.", myN, ".",200,".Rdata"))
mycis <- arrays

# Fagerland figure 5 (left panel)
myN <- 25
# Fagerland book, figure 8.12
myN <- 30
#system.time(mycis <- cifun(n=myN, contrast="RD", alph = c(0.05)))[[3]]/60

#p0 <- as.numeric(dimnames(arrays$mastercp)[[1]])
del <- 0.3
phi <- 0
psi <- 2
#length(p0)
p0 <- seq(0, 1 - del, length.out = 51)
p2 <- p0
#p2[1] <- 0.0001

# Check Newcombe Table 8.5
# - confirm surprisingly shorter widths for MOVER-NW with large phi
myN <- 10
del <- 0
phi <- 0.96
dimnames(arrays$mastercp)
arrays$mastercp["0.4975", "0.4975", "0.96", , "95", "len", ,]
arrays$mastercp["0.4975", "0.4975", "0.96", , "95", "avecp", ,]


p1 <- p2 + del
cp1 <- onecpfun(
  p1 = p1,
  p2 = p2,
  ciarrays = mycis,
  alph = 0.05,
#  phis = phi
  psis = psi
)
widthteam <- c("AS", "SCAS", "SCAS-bc")
widthteam <- c("MOVER-NW", "MOVER-NJ", "BP")
widthteam <- c("AS", "SCAS-bc", "MOVER-NJ", "BP")
widthteam <- c("AS", "BP", "MOVER-NW", "SCAS-bc") # Selected methods from Fagerland plot, plus SCAS-bc
par(pty = "s")
plot(p2,
     cp1[,"AS","len"],
     type = "n",
     ylim = c(0.35, 0.55),
     ylab = "Expected width",
     xlab = "p2",
#     main = paste0("N = ", myN, ", θ = ", del, ", ϕ = ", phi)
     main = paste0("N = ", myN, ", θ = ", del, ", ψ = ", psi)
)
#widthteam <- c("SCAS-bc", "AS", "MOVER-NW", "MOVER-NJ", "BP")
lwds <- c(2, 2, 1, 2, 1, 1)
ltys <- c(1, 2, 0, 3)
mysymbols <- c(NA,NA,5,NA)
for (i in 1:length(widthteam)) {
  lines(x = p2,
        y = cp1[, widthteam[i], "len"],
        lty = ltys[i],
        lwd = lwds[i])
  points(x = p2,
        y = cp1[, widthteam[i], "len"],
        pch = mysymbols[i],
        lwd = lwds[i])
}
legend(x = "bottom", legend = widthteam, lty = ltys, pch = mysymbols, lwd = lwds)

dev.off()

# 2-D interval width plot, RR
myN <- 40
load(file = paste0(outpath, "cis.RR.", myN,".Rdata"))
mycis <- ciarrays

# ?Fagerland figure 6
myN <- 15
#system.time(mycis <- cifun(n = myN, contrast="RR", alph = c(0.05)))[[3]]/60
dimnames(mycis$cis)
mycis$cis[,,"SCAS-bc",,,]
#p0 <- as.numeric(dimnames(arrays$mastercp)[[1]])
theta <- 6
#phi <- 0.1
psi <- 3
#length(p0)
p0 <- seq(0, 1, length.out = 51)
p1 <- p0[c(-1, -51)]
#p2[1] <- 0.0001

#p2 <- p0[p0 + del >= 0 & p0 + del <= 1]

p2 <- p1 / theta
cp1 <- onecpfun(
  p1 = p1,
  p2 = p2,
  ciarrays = mycis,
  alph = 0.05,
#  phis = phi
  psis = psi
)
aslen <- cp1[, "AS" , "len"]
scaslen <- cp1[, "SCAS" , "len"]
(exp(scaslen) - exp(aslen)) / exp(aslen)


#widthteam <- c("AS", "SCAS", "SCAS-bc") #, "MOVER-NW", "MOVER-NJ", "BP")
widthteam <- c("SCAS-bc", "SCAS", "AS", "MOVER-NW", "MOVER-NJ", "BP-W")
widthteam <- c("SCAS-bc", "AS", "MOVER-NW", "MOVER-NJ") #, "BP-W")
#widthteam <- c("AS", "MOVER-NW", "MOVER-NJ", "BP-W")
widthteam <- c("AS", "BP", "MOVER-W", "Wald", "SCAS", "SCAS-bc")
par(pty = "s")
plot(p1,
     cp1[,"AS","len"],
     type = "n",
     ylim = c(2, 6),
     ylab = "Expected width (log)",
     xlab = "p1",
#     main = paste0("N = ", myN, ", θ = ", theta, ", ϕ = ", phi)
     main = paste0("N = ", myN, ", θ = ", theta, ", ψ = ", psi)
)
lwds <- c(2, 2, 1, 2, 1, 1)
ltys <- c(1, 2, 0, 3, 0, 5, 3)
mysymbols <- c(NA,NA,3,NA, 1)
for (i in 1:length(widthteam)) {
  lines(x = p1,
        y = cp1[, widthteam[i], "len"],
        lty = ltys[i],
        lwd = lwds[i])
  points(x = p1,
         y = cp1[, widthteam[i], "len"],
         pch = mysymbols[i],
         lwd = lwds[i])
}
legend(x = "topright", legend = widthteam, lty = ltys, pch = mysymbols, lwd = lwds)
#cp1[, "AS" , "len"]



load(file=paste0(outpath, "cparrays.RD.", 40, ".",200,".Rdata"))
mycps <- arrays$mastercp[, , "0.1", c("AS", "SCAS", "MOVER-W"), "95","cp",,]

dimnames(arrays$mastercp)
#, "SCAS", "MOVER-W"
arrays$mastercp[100:110, 100:110 , "0.1", c("AS", "SCAS", "MOVER-W"), "95","cp",,]
p1 <- p2 <-  as.numeric(dimnames(mycp)[[1]])
del <- 0.2

p1diag <- paste(p1[p2 + del >= 0 & p2 + del <= 1])
p2diag <- paste(p2[p1 - del >= 0 & p1 - del <= 1])
#p1diag <- paste(as.numeric(p1diag) - del)
#diag(mycp[p1diag, p2diag, ]

diagcp <- diag(mycps[p1diag, p2diag, "SCAS"])
length(diagcp)

plot(p2diag,
     diagcp,
     type = "l",
     ylim = c(0.92,1),
     ylab = "CP",
     xlab = "p2",
     main = "CP")
#     main=paste("Coverage probability for NJ,\n at increments of delta from 0 to 0.3. (n1=n2=40)")) #,
abline(h=0.95)
abline(h=0.945,lty=2)

dels <- del+seq(-0.01,0.01,0.005)
for(i in 1:length(dels)){
  deli <- dels[i]
#  p1diag <- paste(p1[round(p1 - deli, 3) >= 0 & round(p1 - deli, 3) <= 1])
#  p2diag <- paste(p2[round(p2 + deli, 3) >= 0 & round(p2 + deli, 3) <= 1])
  p2diag <- paste(p2[p1 - deli >= 0 & p1 - deli <= 1])
  p1diag <- paste(p1[p2 + deli >= 0 & p2 + deli <= 1])
#  p1diag <- p2diag - deli
  diagcpi <- diag(mycps[p1diag, p2diag, "SCAS"])
  lines(p2diag, diagcpi, type="l", lty=2)
}






#### unpaired 2-D plot

#for cross-checking against 2-d plots in papers, e.g. Brown&Li, Fagerland:
root<-"/Users/ssu/Documents/"
#outpath=paste(root,"Main/Courses_papers/skewscore/plots/",sep="")
#outpath=paste(root,"Main/Contract/ssu2010_003/DiffBin/plots/",sep="")
#path=paste(root,"Main/Contract/ssu2010_003/DiffBin/R/diffBinconf/",sep="")
#path2=paste(root,"Main/Contract/ssu2010_003/DiffBin/R/",sep="")
newpath=paste(root,"Main/Courses_papers/skewscore/R/",sep="")
source(paste(newpath,'diffBinconf.all.ssc.R',sep=""))


m=39;n=11
p.grid<-ppoints(100)
p.grid<-seq(0,1,0.01)
#p2<-c(0.2,0.5) #horizontal slice
p2<-p.grid #diagonal slice
p1<-p.grid #add random jitter here if required
methods<-dimnames(diffBinconf.all(5,5,10,10))[[1]]
mastercp=array(NA,dim=c(length(p1),length(p2),length(methods),4,6,1))
dimnames(mastercp)<-list(paste(p1),paste(p2),methods,c(80,90,95,99),c("cp","lncp","avecp","avelncp","len","pow"))
summaries<-array(NA,dim=c(length(methods),5,22,1))
dimnames(summaries)=list(methods,c(80,90,95,98,99),c("meanlen","meanCP","minCP","pctCons","pctBad","pctnear","pctAveCons","pctAvenear","pctAvenearlo","pctAvenearhi","pctAveBad","MNCP","DNCP","maxLNCP","pctBad.1side","pctCons.1side","pctnear.1side","pctAvenear.1side","typeI","maxtypeI","avetypeI","maincolour"))
cis<-array(NA,dim=c(length(methods),2,(m+1)*(n+1),4,1))
dimnames(cis)[[4]]<-c(80,90,95,99)
cis<-array(NA,dim=c(length(methods),3,(m+1)*(n+1),5,1))
dimnames(cis)[[4]]<-c(80,90,95,98,99)
dimnames(cis)[[2]]<-c("LCL","UCL","disjoint")
dim(cis)
myarrays<-list(mastercp,summaries,cis)
rm(mastercp)
#rm(mastercp.Unif)
rm(summaries)
system.time(slice<-bigfun(m=m,n=n,mm=m,nn=n,p1=p1,p2=p2,alph=0.05,arrays=myarrays,square=TRUE,jitt=FALSE))[3]/60




#Fagerland plot & why their conclusion on MN is wrong.  They've used Mee for starters.
del <- 0
fmethod<-"Mee"
attrib<-"cp"
#attrib<-"len"
#attrib<-"pow"
#p1diag <- paste(p1[p2+del>=0 & p2+del<=1])
#p2diag <- paste(p2[p1-del>=0 & p1-del<=1])
p1diag <- paste(p1[p1-del>=0 & p1-del<=1])
p2diag <- paste(p2[p2+del>=0 & p2+del<=1])
diagplot <- diag(slice[[1]][p1diag,p2diag,fmethod,"95",attrib,1])
diagdata <- (slice[[1]][p1diag,p2diag,,"95",attrib,1])
length(diagplot)
#plot(p2diag,diagplot,type="l",ylim=c(0.92,1),ylab="power",xlab="p2",main=paste("Coverage probability for NJ,\n at increments of delta from 0 to 0.3. (n1=n2=40)")) #,
#plot(p2diag,diagplot,type="l",ylim=c(0,0.05),ylab="power",xlab="p2",main=paste("Coverage probability for NJ,\n at increments of delta from 0 to 0.3. (n1=n2=40)")) #,
#plot(p1diag,diagplot,type="l",ylim=c(0.95,0.98),ylab="power") #,

plot(p2diag,diagplot,type="l",ylim=c(0.92,1),ylab="CP",xlab="p2",
     main=paste0("Coverage probability for Mee/MN/SCAS,\n at delta=", del, " (n1= ", m, " n2=", n, ")")) #,

abline(h=0.95)
abline(h=0.945,lty=2)

#main=expression(paste("Power of Brown-Li CI for p1=p2, m=n=20")))
lines(p2diag,diag(diagdata[,,"MN"]),lty=1)
lines(p2diag,diag(diagdata[,,"Mee"]),lty=2, lwd=2)
lines(p2diag,diag(diagdata[,,"SC"]),lty=3, lwd=4)
#lines(p2diag,diag(diagdata[,,"N"]),lty=4, lwd=1)
#lines(p2diag,diag(diagdata[,,"Wald"]),lty=4)






# From the archive

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


# For troubleshooting plots

load(file = paste0(outpath, "cparrays.RD.", 20, ".",200,".Rdata"))
CPcontour(
  plotdata = arrays,
                      alpha = 0.05,
                      par3 = 0.75,
                      nums = "20",
                      xlim = c(0,1),
                      ylim = c(0,1),
                      methlab = "MOVER-NJ",
                      avg = T,
                      lside = F,
                      lines = F,
                      lines1 = F,
                      CIlen = T,
                      locind = F,
                      crude = F,
                      res.factor = 6,
                      colour = T,
                      textsize = 1

                      )


}
