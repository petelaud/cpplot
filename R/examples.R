if (FALSE) {
x <- c(1, 1, 7, 12)

# TABLE 5.1 Example CIs for RD
round(rbind(
  pairbinci(x = x, contrast = "RD", method_RD = "Score", bcf = TRUE, skew = TRUE)$estimates[,1:3], # SCAS-bc
  pairbinci(x = x, contrast = "RD", method_RD = "Score", bcf = FALSE, skew = TRUE)$estimates[,1:3], # SCAS(no bc)
  pairbinci(x = x, contrast = "RD", method_RD = "Score", bcf = FALSE, skew = FALSE)$estimates[,1:3], # AS
  pairbinci(x = x, contrast = "RD", method_RD = "MOVER_newc", moverbase = 'jeff')$estimates, # MOVER-NJ
  pairbinci(x = x, contrast = "RD", method_RD = "MOVER_newc", moverbase = 'wilson')$estimates, # MOVER-W
  pairbinci(x = x, contrast = "RD", method_RD = "MOVER", moverbase = 'wilson')$estimates, # MOVER-W
  pairbinci(x = x, contrast = "RD", method_RD = "BP")$estimates, # BP
  pairbinci(x = x, contrast = "RD", method_RD = "Score", bcf = TRUE, skew = TRUE, cc=0.125)$estimates[,1:3], # SCAS-cc
  pairbinci(x = x, contrast = "RD", method_RD = "Score", bcf = TRUE, skew = TRUE, cc=0.5)$estimates[,1:3], # SCAS-cc
  pairbinci(x = x, contrast = "RD", method_RD = "MOVER_newc", moverbase = 'jeff', cc=0.125)$estimates, # MOVER-NJcc
  pairbinci(x = x, contrast = "RD", method_RD = "MOVER_newc", moverbase = 'jeff', cc=0.5)$estimates # MOVER-NJcc
), 3)[, c(1, 3)]

# TABLE 5.2 Example CIs for RR
round(rbind(
  pairbinci(x = x, contrast = "RR", method_RR = "Score", bcf = TRUE, skew = TRUE)$estimates[,1:3], # SCAS-bc
  pairbinci(x = x, contrast = "RR", method_RR = "Score", bcf = FALSE, skew = TRUE)$estimates[,1:3], # SCAS(no bc)
  pairbinci(x = x, contrast = "RR", method_RR = "Score", bcf = FALSE, skew = FALSE)$estimates[,1:3], # AS
  pairbinci(x = x, contrast = "RR", method_RR = "MOVER_newc", moverbase = 'jeff')$estimates, # MOVER-NJ
  pairbinci(x = x, contrast = "RR", method_RR = "MOVER_newc", moverbase = 'wilson')$estimates, # MOVER-W
  pairbinci(x = x, contrast = "RR", method_RR = "MOVER", moverbase = 'wilson')$estimates, # MOVER-W
  pairbinci(x = x, contrast = "RR", method_RR = "BP", moverbase = 'wilson')$estimates, # BP
  pairbinci(x = x, contrast = "RR", method_RR = "Score", bcf = TRUE, skew = TRUE, cc=0.125)$estimates[,1:3], # SCAS-cc
  pairbinci(x = x, contrast = "RR", method_RR = "Score", bcf = TRUE, skew = TRUE, cc=0.5)$estimates[,1:3], # SCAS-cc
  pairbinci(x = x, contrast = "RR", method_RR = "MOVER_newc", moverbase = 'jeff', cc=0.125)$estimates, # MOVER-NJcc
  pairbinci(x = x, contrast = "RR", method_RR = "MOVER_newc", moverbase = 'jeff', cc=0.5)$estimates # MOVER-NJcc
), 3)[, c(1, 3)]



# Equivalence illustration

x <- c(240, 37, 18, 40)
x <- c(240, 30, 11, 40)
x <- c(240, 31, 12, 40)
x <- c(217, 30, 13, 40)
sum(x)
pairbinci(x = x, contrast = "RD", method_RD = "Score", bcf = TRUE, skew = TRUE)$estimates[,1:3] # SCAS-bc
pairbinci(x = x, contrast = "RD", method_RD = "Score", bcf = TRUE, skew = FALSE)$estimates[,1:3] # SCAS-bc
pairbinci(x = x, contrast = "RD", method_RD = "MOVER")$estimates[,1:3] #




}
