#======================#
#                      #
##### DUKIC METHOD #####
#                      #
#======================#

#############################
### SCRIPT SETUP COMMANDS ###
#############################

### Clear Everything ###
rm(list = ls())

### Load Required Libraries ###
library(MASS)

### Set Directory ###
data.dir <- "~/Data/Kathy" # Where the data is.
dukic.dir <- "~/Documents/Git/Dukic-Method" # Where the dukic-method.R script is.
output.dir <- "~/Desktop" # Where you want your graphs output.

#####################################
### SOURCE DUKIC-METHOD FUNCTIONS ###
#####################################

setwd(dukic.dir)
source("dukic-method.R")

#######################################
### CREATE DATASETS & DATA CLEANING ###
#######################################

### Read Data Into R ###
# The programme expects a data frame called "data" with the following 8 columns:
# "study", "cut.off", "n", "TP", "FP", "TN", "FN"

setwd(data.dir)
data <- read.csv("newROC1.csv", stringsAsFactors = FALSE)
data <- data[which(names(data) %in% c("study", "cut.off", "n", "TP", "FP", "FN", "TN"))]
cat("Warning: have to replace", sum(data$n!=data$TP+data$FP+data$FN+data$TN), "values of n \n")
data$n <- data$TP + data$FP + data$FN + data$TN # Integriy check.
data$sens <- data$TP / (data$TP + data$FN)
data$spec <- data$TN / (data$TN + data$FP)

#########################################
### GENERATE DIAGNOSTIC ACCURACY PLOT ###
#########################################

### Plot Parameters ###
dukic <- fDukic(data)

solution <- optim(par = c(0.75, 200, 5),
                  fLoglikelihood,
                  method = "Nelder-Mead",
                  control = list(fnscale = -1,
                                 maxit = 1e5,
                                 reltol = 1e-12,
                                 trace = 1L), # Specify maximisation, as opposed to minimisation, through fnscale parameter
                  hessian = TRUE)

if(solution$convergence == 0) {
  print(paste("alpha =", round(solution$par[1], 2)), quote = F)
  print(paste("beta =", round(solution$par[2], 2)), quote = F)
  print(paste("intercept =", round(solution$par[3], 2)), quote = F)
  print(paste("AUC =", round(fAUC(function(x) {fSpectosens(1 - x, intercept = solution$par[3], alpha = solution$par[1], beta = solution$par[2])}, a = 0, b = 1, verbose = F), 3)), quote = F)
}

x <- seq(1 - max(data$spec), 1 - min(data$spec), 0.001)
y <- fSpectosens(1 - x, solution$par[3], solution$par[1], solution$par[2])

### Plot SROC ###
setwd(output.dir)
pdf("plots.pdf")

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = 1)
axis(side = 1, at = seq(0, 1, 0.1))
axis(side = 2, at = seq(0, 1, 0.1))
title(xlab = "False Positive Rate (1 - Specificity)", ylab = "True Positive Rate (Sensitivity)")
lines(y ~ x, type = "l", col="black")

### Annotate Plot ###
sens <- data$TP / (data$TP + data$FN)
spec <- 1 - data$FP / (data$FP + data$TN)
spec90 <- fSenstospec(0.9, solution$par[3], solution$par[1], solution$par[2])
cat("BNP: when sensitivity == 0.9, specificity == ", spec90, "\n")
cat("LR- at this point is ", 0.1 / spec90, '\n')
sens90 <- fSpectosens(0.9, solution$par[3], solution$par[1], solution$par[2])
cat("BNP: when specificity == 0.9, sensitivity == ", sens90, "\n")
cat("LR+ at this point is ", sens90 / 0.1, '\n')
points(1 - spec, sens, pch = 16, bg = "black", col = "black", cex = 0.5)
#text(1 - spec, sens + 0.015, labels = data$cut.off, col = "black", cex = 0.5)
#text(1 - spec, sens + 0.03, labels = data$study, col = "black", cex = 0.25)

### Generate Confidence Intervals for Alpha ###
solution$var <- solve(-solution$hessian)
print(cbind( solution$par, solution$par - 1.96 * sqrt(diag(solution$var)), solution$par + 1.96 * sqrt(diag(solution$var)) ))

### Plot SROC Confidence Intervals ###
cat("Theta", "\t", "Sensitivity (CI) [Range]", "\t", "Specifity (CI) [Range]", "\n")
xoffset <- c(0, 0, 0)	# xoffset and ymult just tweak the position of the numerical labels
#points <- quantile(data$cut.off, seq(0.25, 0.5, 0.75))
points <- c(100, 200, 300)
for (i in 1:length(points)) {
  sens <- fThetatosens(points[i], solution$par[3], solution$par[1], solution$par[2])
  spec <- fThetatospec(points[i], solution$par[3], solution$par[1], solution$par[2])
  points(1 - spec, sens, col = "black", pch = 3)
  lSim <- fSimulateCI(points[i], solution$par, solution$hessian)
  fPlotsimci(lSim, "black")
  text(1 - lSim$lSpec$mle + xoffset[i], lSim$lSens$upper + 0.02, labels = paste("BNP<", points[i], sep = ""), col = "black", cex = 0.7)
  cat(formatC(lSim$theta, digits = 2, format = "f"), '\t', sep = "", file = "")
  fTabulate(lSim$lSens, data[which(data$cut.off == points[i]), ]$sens)
  cat("\t")
  fTabulate(lSim$lSpec, data[which(data$cut.off == points[i]), ]$spec)
  cat("\n")
}

cat("\n")
cat("\nTheta", "\t", "LR+ (CI) [Range]", "\t\t\t\t\t", "LR- (CI) [Range]", "\n")
data$lrplus <- data$sens / (1 - data$spec)
data$lrminus <- (1 - data$sens) / data$spec
for (i in 1:length(points)) {
  lSim <- fSimulateCI(points[i], solution$par, solution$hessian)
  cat(formatC(lSim$theta, digits = 2, format = "f"), '\t', sep = "", file = "")
  fTabulate(lSim$lLrplus, data[which(data$cut.off == points[i]), ]$lrplus)
  cat("\t")
  fTabulate(lSim$lLrminus, data[which(data$cut.off == points[i]), ]$lrminus)
  cat("\n")
}

### ggplot 2 Plot ###

dev.off()