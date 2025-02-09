if (FALSE) {
# Plot of type I error for SCAS test ('N-1' adjusted McNemar test)

  # AS test matches McNemar
  pairbinci(x = c(1,1,7,12), bcf = FALSE)$pval[, 2]
  mcnemar.test(x = matrix(c(1,1,7,12), nrow=2), correct=FALSE)
  # 'N-1' AS test:
  pairbinci(x = c(1,1,7,12), bcf = TRUE)$pval[, 2]
  # conditional mid-p test
  midp <- 2*pbinom(1, 8, 0.5, lower.tail = TRUE) - dbinom(1, 8, 0.5)

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
#     ylab = "Type I error rate",
#     xlab = "p2 = p1",
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
#lines(p2, 1 - cp1[,"TDAS","cp"], lty=3, lwd = 2)
#lines(p2, 1 - cp1[,"SCASstrat","cp"], lty=4, lwd = 2)
dev.off()






# Power
p1 <- 0.35
p2 <- 0.65
psi <- 2
p1 <- 0.1
p2 <- 0.4
psi <- 2
n <- 38
g <- (expand.grid(x11 = 0:n, x12 = 0:n, x21 = 0:n))
# reduce to possible combinations of a,b,c for paired data.
g <- g[(g$x12 <= n - g$x11) &
         (g$x21 <= n - g$x11 - g$x12), ]
xs <- data.matrix(cbind(g, x22 = n - rowSums(g)))
prob <- pdfpair(p1 = p1,
                p2 = p2,
                psi = psi,
                x = xs)
pvals <- sapply(1:dim(xs)[[1]], function(i) pairbinci(xs[i,], contrast="RD", method_RD="Score_closed", skew=FALSE, bcf=FALSE)$pval[2])
powercalc <- (pvals < 0.05) %*% prob
powercalc

pvals <- sapply(1:dim(xs)[[1]], function(i) pairbinci(xs[i,], contrast="RD", method_RD="Score", skew=TRUE, bcf=TRUE)$pval[2])
powercalc <- (pvals < 0.05) %*% prob
powercalc


# Repeat for the following:
#Power was calculated for N = 1, 2, . . . , 100, θ = 1.0, 2.0,
#3.0, 5.0, 10.0, p1+ = 0.1, 0.35, 0.6, and   = p+1 − p1+ =
#  0.10, 0.15, 0.20, 0.25, 0.30, 0.35.

# Reproduce Fagerland2013 evaluation of Type I error rate (TIER)

n <- 10
psi <- c(2)
p1 <- p2 <- 0.1
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
#  sum(prob)
  xsub <- xs[prob > 1E-8, ]
#  dim(xsub)
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
#    tier <- (pvals == pvals) %*% prob[prob > 1E-8]
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
  c(tier, tier2)
}

myparams <- expand.grid(p1 = seq(0, 1, 0.01),
                        psi = c(1, 2, 3, 5, 10),
                        n = seq(10, 100, 5)
)
myparams <- expand.grid(p1 = seq(0, 1, 0.02),
                        phi = seq(0.25, 0.75, 0.05),
                        n = seq(10, 100, 5)
)
dim(myparams)

# Runtime: 3.3 hours
system.time(tiers <- sapply(1:dim(myparams)[1], function(i) tier(rev(unlist(myparams[i,])))))[[3]]/60

system.time(tiers <- sapply(1:dim(myparams)[1], function(i) tier(rev(unlist(myparams[i,])))))[[3]]/60
#tiers
apply(tiers, 1, mean)
apply(tiers, 1, max)
apply(tiers, 1, function(x) mean(x > 0.05))
apply(tiers, 1, function(x) mean(x < 0.03))
mytiers <- t(tiers)
dim(mytiers)
dimnames(mytiers)[[2]] <- c("SCAStier", "AStier")
mytiers <- cbind(myparams, mytiers)
save(mytiers, file = paste0(outpath, "tiers.Rdata"))


par(mfrow=c(1,1), pty='s', mar=(c(2,4,2,0)+0.1))

sel <- (myparams[,2] == 1 & myparams[, 3] == 40)
subparams <- myparams[sel, ]
plot(
  subparams[, 1], tiers[1, sel], type='n', ylim=c(0, 0.06))
abline(h=0.05)
lines(
  subparams[, 1], tiers[1, sel])
lines(
  subparams[, 1], tiers[2, sel], lty=2)

mean(tiers[1, sel])

#scorepair()

}
