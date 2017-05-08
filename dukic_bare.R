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
data.dir <- "~/Downloads" # Where the data is.
dukic.dir <- "~/Downloads" # Where the dukic-method.R script is.
output.dir <- "~/Downloads" # Where you want your graphs output.

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
data <- read.csv("ROC1.csv", stringsAsFactors = FALSE)
data <- data[which(names(data) %in% c("study", "cut.off", "n", "TP", "FP", "FN", "TN"))]
cat("Warning: have to replace", sum(data$n!=data$TP+data$FP+data$FN+data$TN), "values of n \n")
data$n <- data$TP + data$FP + data$FN + data$TN # Integriy check.

#########################################
### GENERATE DIAGNOSTIC ACCURACY PLOT ###
#########################################

### Plot Parameters ###
x <- c(0:100) / 100

dukic <- fDukic(data)
solution <- optim(par = c(0, 1, 20),
                  fLoglikelihood,
                  control = list(fnscale = -1), # Specify maximisation, as opposed to minimisation, through fnscale parameter
                  hessian = TRUE)
print(solution$par)

### Plot SROC ###
setwd(output.dir)
pdf("plots.pdf")

y <- fSpectosens(1 - x, solution$par[3], solution$par[1], solution$par[2])
plot(y ~ x, type = "l", col="black", xlab = "1-Specificity", ylab = "Sensitivity", xaxs = 'i', yaxs = 'i')

### Annotate Plot ###
sens <- data$TP / (data$TP + data$FN)
spec <- 1 - data$FP / (data$FP + data$TN)
spec90 <- fSenstospec(0.9, solution$par[3], solution$par[1], solution$par[2])
cat("PCT: when sensitivity == 0.9, specificity == ", spec90, "\n")
cat("LR- at this point is ", 0.1 / spec90, '\n')
sens90 <- fSpectosens(0.9, solution$par[3], solution$par[1], solution$par[2])
cat("PCT: when specificity == 0.9, sensitivity == ", sens90, "\n")
cat("LR+ at this point is ", sens90 / 0.1, '\n')
points(1 - spec, sens, pch = 16, bg = "black", col = "black", cex = 1.3)
text(1 - spec, sens + 0.015, labels = data$cut.off, col = "black", cex = 0.67)
text(1 - spec, sens + 0.03, labels = data$study, col = "black", cex = 0.5)

### Generate Confidence Intervals for Alpha ###
solution$var <- solve(-solution$hessian)
print(cbind( solution$par, solution$par - 1.96 * sqrt(diag(solution$var)), solution$par + 1.96 * sqrt(diag(solution$var)) ))

### Plot SROC Confidence Intervals ###
cat("Theta", "\t", "Sensitivity (CI) [Range]", "\t", "Specifity (CI) [Range]", "\n")
data$sens <- data$TP / (data$TP + data$FN)
data$spec <- data$TN / (data$TN + data$FP)
xoffset <- c(0, 0, 0)	# xoffset and ymult just tweak the position of the numerical labels
points <- sort(unique(data$cut.off))
for (i in 1:length(points)) {
  sens <- fThetatosens(points[i], solution$par[3], solution$par[1], solution$par[2])
  spec <- fThetatospec(points[i], solution$par[3], solution$par[1], solution$par[2])
  points(1 - spec, sens, col = "black", pch = 3)
  lSim <- fSimulateCI(points[i], solution$par, solution$hessian)
  fPlotsimci(lSim, "black")
  text(1 - lSim$lSpec$mle + xoffset[i], lSim$lSens$upper + 0.02, labels = paste("PCT<", points[i], sep = ""), col = "black", cex = 0.7)
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
dev.off()