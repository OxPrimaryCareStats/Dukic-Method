#########################################################################
# ===================== Van den Bruel et al, 2011 ===================== #
#########################################################################

# ===================== Read in Data Frame ===================== #

# Clear memory and load required libraries
rm(list=ls())
require(MASS) # Run the command "install.packages('MASS')" if this package is not installed.

# Read data frame into memory
data<-read.csv('vandenbruel.csv')
cat("Warning: have to replace", sum(data$n!=data$TP+data$FP+data$FN+data$TN), "values of n \n")
data$n<-data$TP+data$FP+data$FN+data$TN

# Where to show CIs
#cipoints<-c(0.6, 0.15, 0.03)
#cioffset<-0.003

# ===================== Load Dukic Method Functions ===================== #

source('dukic-method.R')

# ===================== Plot SROC Values ===================== #

# Plot parameters
x<-c(0:100)/100

print(" --------------------- PCT ------------------------------------------------")
pctdata<-data[which(data$test=="PCT"),]
dukic<-pctdukic<-fDukic(pctdata)
solution<-optim(par=c(0,0,0), fLoglikelihood, control=list(fnscale=-1), hessian=TRUE) # Specify maximisation, as opposed to minimisation, through fnscale parameter
print("Solution for PCT:")
print(solution$par)

# Set up plotting characters
initial<-substr(pctdata$study,1,1)
pctpch<-21*(initial=="A") + 22*(initial=="B") + 23*(initial=="G") + 24*(initial=="H") + 25*(initial=="T")

# Plot SROC
y<-fSpectosens(1-x, solution$par[3], solution$par[1], solution$par[2])
plot(y~x, type="l", col="black", xlab="1-Specificity", ylab="Sensitivity", xaxs='i', yaxs='i')
sens<-pctdata$TP/(pctdata$TP+pctdata$FN)
spec<-1-pctdata$FP/(pctdata$FP+pctdata$TN)
pctsolution<-solution	# Save results from PCT analysis
spec90<-fSenstospec(0.9, solution$par[3], solution$par[1], solution$par[2])
cat("PCT: when sensitivity==0.9, specificity==", spec90, "\n")
cat("LR- at this point is ", 0.1/spec90, '\n')
# print(fLRminusCI(0.9, solution$par, solution$hessian))
sens90<-fSpectosens(0.9, solution$par[3], solution$par[1], solution$par[2])
cat("PCT: when specificity==0.9, sensitivity==", sens90, "\n")
cat("LR+ at this point is ", sens90/0.1, '\n')
# print(fLRplusCI(0.9, solution$par, solution$hessian))
points(1-spec, sens, pch=16, bg="black", col="black", cex=1.3)
text(1-spec, sens+0.015, labels=pctdata$cut.off, col="black", cex=0.67)
text(1-spec, sens+0.03, labels=pctdata$study, col="black", cex=0.5)

print("--------------------- CRP ------------------------------------------------")
# Repeat for CRP
crpdata<-data[which(data$test=="CRP"),]
# Set up plotting characters
initial<-substr(crpdata$study,1,1)
crppch<-21*(initial=="A") + 22*(initial=="B") + 23*(initial=="G") + 24*(initial=="H") + 25*(initial=="T")
sens<-crpdata$TP/(crpdata$TP+crpdata$FN)
spec<-1-crpdata$FP/(crpdata$FP+crpdata$TN)
sens<-pmin(sens, 0.99)	# Over-writing the sens=100% point
points(1-spec, sens, pch=1, bg="red", col="red", xlim=c(0,1), ylim=c(0,1), xlab="1-specificity", ylab="sensitivity",  cex=1.3)
text(1-spec, sens-0.015, labels=crpdata$cut.off, col="red", cex=0.67)
text(1-spec, sens-0.03, labels=crpdata$study, col="red", cex=0.5)

dukic<-crpdukic<-fDukic(crpdata)
solution<-optim(par=c(1,1,1), fLoglikelihood, control=list(fnscale=-1), hessian=TRUE) # Specify maximisation, as opposed to minimisation, through fnscale parameter
y<-fSpectosens(1-x, solution$par[3], solution$par[1], solution$par[2])
lines(x, y, type="l", col="red", lty="dotted")

crpsolution<-solution	# Save results from CRP analysis
spec90<-fSenstospec(0.9, solution$par[3], solution$par[1], solution$par[2])
cat("CRP: when sensitivity==0.9, specificity==", spec90, "\n")
cat("LR- at this point is ", 0.1/spec90, '\n')
sens90<-fSpectosens(0.9, solution$par[3], solution$par[1], solution$par[2])
cat("CRP: when specificity==0.9, sensitivity==", sens90, "\n")
cat("LR+ at this point is ", sens90/0.1, '\n')

print("======== Confidence intervals for alpha ===========")
print("-------- PCT -----------")
pctsolution$var<-solve(-pctsolution$hessian)
print(cbind( pctsolution$par, pctsolution$par-1.96*sqrt(diag(pctsolution$var)), pctsolution$par+1.96*sqrt(diag(pctsolution$var)) ))
print("-------- CRP -----------")
crpsolution$var<-solve(-crpsolution$hessian)
print(cbind( crpsolution$par, crpsolution$par-1.96*sqrt(diag(crpsolution$var)), crpsolution$par+1.96*sqrt(diag(crpsolution$var)) ))

# legend
legend(x=0.5,y=0.3, c("PCT data ", "CRP data ", "PCT estimate", "CRP estimate", "95% Confidence intervals"), col=c("black", "red", "black", "red", "black"), pch=c(16,1, 3, 3, 3), lty=c("blank", "blank", "solid", "dotted", "dashed"))

# ===================== Plot SROC Confidence Intervals ===================== #

cat("Theta", "\t", "Sensitivity (CI) [Range]", "\t", "Specifity (CI) [Range]", "\n")
print(" --------------------- PCT ------------------------------------------------")
pctdata$sens<-pctdata$TP / (pctdata$TP + pctdata$FN)
pctdata$spec<-pctdata$TN / (pctdata$TN + pctdata$FP)
xoffset<-c(0,0,0)	# xoffset and ymult just tweak the position of the numerical labels
pctpoints<-c(0.5, 1, 2)
for (i in 1:length(pctpoints)) {
  sens<-fThetatosens(pctpoints[i], pctsolution$par[3], pctsolution$par[1], pctsolution$par[2])
  spec<-fThetatospec(pctpoints[i], pctsolution$par[3], pctsolution$par[1], pctsolution$par[2])
  points(1-spec, sens, col="black", pch=3)
  lSim<-fSimulateCI(pctpoints[i], pctsolution$par, pctsolution$hessian)
  fPlotsimci(lSim, "black")
  text(1-lSim$lSpec$mle+xoffset[i], lSim$lSens$upper+0.02, labels=paste("PCT<", pctpoints[i], sep=""), col="black", cex=0.7)
  cat(formatC(lSim$theta, digits=2, format="f"), '\t', sep="", file="")
  fTabulate(lSim$lSens, pctdata[which(pctdata$cut.off==pctpoints[i]),]$sens)
  cat("\t")
  fTabulate(lSim$lSpec, pctdata[which(pctdata$cut.off==pctpoints[i]),]$spec)
  cat("\n")
}

print("--------------------- CRP ------------------------------------------------")
crpdata$sens<-crpdata$TP / (crpdata$TP + crpdata$FN)
crpdata$spec<-crpdata$TN / (crpdata$TN + crpdata$FP)
xoffset<-c(-0.02,-0.02,0.005)	# xoffset and ymult just tweak the position of the numerical labels
ymult<-c(1,1,-1)
crppoints<-c(20, 40, 80 )
for (i in 1:length(crppoints)) {
  sens<-fThetatosens(crppoints[i], crpsolution$par[3], crpsolution$par[1], crpsolution$par[2])
  spec<-fThetatospec(crppoints[i], crpsolution$par[3], crpsolution$par[1], crpsolution$par[2])
  points(1-spec, sens, col="red", pch=3)
  lSim<-fSimulateCI(crppoints[i], crpsolution$par, crpsolution$hessian)
  fPlotsimci(lSim, "red")
  text(1-lSim$lSpec$lower+xoffset[i], lSim$lSens$mle+(lSim$lSens$mle-lSim$lSens$lower+0.02)*ymult[i], labels=paste("CRP<",crppoints[i], sep=""), col="red", cex=0.7)
  cat(formatC(lSim$theta, digits=2, format="f"), '\t', sep="", file="")
  fTabulate(lSim$lSens, crpdata[which(crpdata$cut.off==crppoints[i]),]$sens)
  cat("\t")
  fTabulate(lSim$lSpec, crpdata[which(crpdata$cut.off==crppoints[i]),]$spec)
  cat("\n")
}

cat("\n")
cat("\nTheta", "\t", "LR+ (CI) [Range]", "\t\t\t\t\t", "LR- (CI) [Range]", "\n")
print(" --------------------- PCT ------------------------------------------------")
pctdata$lrplus<-pctdata$sens/(1-pctdata$spec)
pctdata$lrminus<-(1-pctdata$sens)/pctdata$spec
for (i in 1:length(pctpoints)) {
  lSim<-fSimulateCI(pctpoints[i], pctsolution$par, pctsolution$hessian)
  cat(formatC(lSim$theta, digits=2, format="f"), '\t', sep="", file="")
  fTabulate(lSim$lLrplus, pctdata[which(pctdata$cut.off==pctpoints[i]),]$lrplus)
  cat("\t")
  fTabulate(lSim$lLrminus, pctdata[which(pctdata$cut.off==pctpoints[i]),]$lrminus)
  cat("\n")
}

print("--------------------- CRP ------------------------------------------------")
crpdata$lrplus<-crpdata$sens/(1-crpdata$spec)
crpdata$lrminus<-(1-crpdata$sens)/crpdata$spec
for (i in 1:length(crppoints)) {
  lSim<-fSimulateCI(crppoints[i], crpsolution$par, crpsolution$hessian)
  cat(formatC(lSim$theta, digits=2, format="f"), '\t', sep="", file="")
  fTabulate(lSim$lLrplus, crpdata[which(crpdata$cut.off==crppoints[i]),]$lrplus)
  cat("\t")
  fTabulate(lSim$lLrminus, crpdata[which(crpdata$cut.off==crppoints[i]),]$lrminus)
  cat("\n")
}
