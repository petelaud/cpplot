# cpslice function to reproduce 2-D plots e.g. in Fagerland et al.

if (FALSE) {



load(file=paste0(outpath, "cparrays.RD.", 65, ".",200,".Rdata"))
p0 <- as.numeric(dimnames(arrays$mastercp)[[1]])
del <- 0.2
del <- 0
length(p0)
p0 <- seq(0,1-del,length.out=51)
p2 <- p0

#p2 <- p0[p0 + del >= 0 & p0 + del <= 1]

# 2-D coverage plot
p1 <- p2 + del
  cp1 <- onecpfun(
    p1 = p1,
    p2 = p2,
    ciarrays = arrays,
    alph = 0.05,
    phis = 0.25
  )
plot(p2,
     cp1[,"AS","cp"],
     type = "l",
     ylim = c(0.90,1),
     ylab = "CP",
     xlab = "p2",
     main = "CP"
)
abline(h=0.95)
rect(
  xleft = par("usr")[1], xright = par("usr")[2], ybottom = 0.94, ytop = 0.95,
  border = NA, col = adjustcolor("gray", alpha = 0.3)
)
#abline(h=0.945,lty=2)

dels <- del + seq(-0.003,0.003,0.001)[-5]
for(i in 1:length(dels)){
  p1 <- p2 + dels[i]
  cp1 <- onecpfun(
    p1 = p1,
    p2 = p2,
    ciarrays = arrays,
    alph = 0.05,
    phis = 0.25
  )
  lines(p2,
       cp1[,"AS","cp"],
       lty = 2
  )
}


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
load(file=paste0(outpath, "cparrays.RD.", 40, ".",200,".Rdata"))
system.time(mycis <- cifun(n=25, contrast="RD", alph = c(0.05)))[[3]]/60
#p0 <- as.numeric(dimnames(arrays$mastercp)[[1]])
del <- 0.3
#length(p0)
p0 <- seq(0, 1 - del, length.out = 51)
p2 <- p0
#p2[1] <- 0.0001

#p2 <- p0[p0 + del >= 0 & p0 + del <= 1]

p1 <- p2 + del
cp1 <- onecpfun(
  p1 = p1,
  p2 = p2,
  ciarrays = mycis,
  alph = 0.05,
  phis = 0.1
)
par(pty = "s")
plot(p2,
     cp1[,"AS","len"],
     type = "n",
     ylim = c(0.3,0.55),
     ylab = "Expected width",
     xlab = "p2",
     main = paste0("N = 25, θ = 0.3, ϕ = 0.1")
)
#widthteam <- c("AS", "SCAS", "SCAS-bc") #, "MOVER-NW", "MOVER-NJ", "BP")
widthteam <- c("SCAS-bc", "AS", "MOVER-NW", "MOVER-NJ", "BP")
lwds <- c(2, 2, 1, 1, 1)
for (i in 1:length(widthteam)) {
  lines(x = p2,
        y = cp1[, widthteam[i], "len"],
        lty = i,
        lwd = lwds[i])
}
legend(x = "bottom", legend = widthteam, lty = 1:length(widthteam), lwd = lwds)



# 2-D interval width plot, RR
load(file=paste0(outpath, "cparrays.RR.", 30, ".",200,".Rdata"))
#system.time(mycis <- cifun(n=30, contrast="RR", alph = c(0.05)))[[3]]/60
#p0 <- as.numeric(dimnames(arrays$mastercp)[[1]])
theta <- 3
#length(p0)
p0 <- seq(0, 1, length.out = 51)
p1 <- p0[c(-1, -51)]
#p2[1] <- 0.0001

#p2 <- p0[p0 + del >= 0 & p0 + del <= 1]

p2 <- p1 / theta
cp1 <- onecpfun(
  p1 = p1,
  p2 = p2,
  ciarrays = arrays,
  alph = 0.05,
  phis = 0.1
)
par(pty = "s")
plot(p1,
     cp1[,"AS","len"],
     type = "n",
     ylim = c(-2, 0),
     ylab = "Expected width",
     xlab = "p1",
     main = paste0("N = 30, θ = 6, ϕ = 0.25")
)
#widthteam <- c("AS", "SCAS", "SCAS-bc") #, "MOVER-NW", "MOVER-NJ", "BP")
widthteam <- c("SCAS-bc", "SCAS", "AS", "MOVER-NW", "MOVER-NJ", "BP-W")
#widthteam <- c("AS", "MOVER-NW", "MOVER-NJ", "BP-W")
lwds <- c(2, 2, 2, 1, 1, 1)
for (i in 1:length(widthteam)) {
  lines(x = p1,
        y = cp1[, widthteam[i], "len"],
        lty = i,
        lwd = lwds[i])
}
legend(x = "topright", legend = widthteam, lty = 1:length(widthteam), lwd = lwds)
cp1[, , "len"]



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


}
