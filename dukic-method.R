############################################################
# ===================== Dukic Method ===================== #
############################################################

# ===================== Read in Data Frame ===================== #

# Read data frame into memory
#data <- read.csv('file_containing__study_data.csv')
# The programme expects a data frame called "data" with the following 8 columns:
# "study", "test", "cut.off", "n", "TP", "FP", "TN", "FN"

# Where to show CIs
#cipoints<-c(0.6, 0.15, 0.03)
#cioffset<-0.003

# ===================== Derivation Functions ===================== #

fLogit <- function(x)
  { log(x / (1 - x)) }

fInvlogit <- function(x)
  { 1 / (1 + exp(-x)) }
stopifnot(fInvlogit(fLogit(0.371)) == 0.371)

fCategorise <- function(n, vTp, vFp, vFn, vTn)
  # Given vectors of true positive tests/false positive tests/etc at each threshold, deduce the number who fell into each test category
  {
  stopifnot( max(length(vTp), length(vFp), length(vFn), length(vTn)) == min(length(vTp), length(vFp), length(vFn), length(vTn)) )
  stopifnot(vTp + vFp + vFn + vTn == n)
  
  k <- length(vTp)	# Length of input vectors; one less than actual # test categories
  Y <- rep(NA, k + 1)	# Will be number (regardless of true status) in each of k+1 categories of test result
  D <- rep(NA, k + 1)	# Will be disease status in each of k+1 test categories
  Y[1] <- vFn[1] + vTn[1]
  D[1] <- vFn[1]
  stopifnot(k >= 2)
  for (j in 2:k){
    Y[j] <- (vFn[j] + vTn[j]) - sum(Y[1:(j - 1)])
    D[j] <- vFn[j] - sum(D[1:(j - 1)])
    }
  Y[k + 1] <- vTp[k] + vFp[k]
  D[k + 1] <- vTp[k]
  if (sum(Y) != n) cat("Warning: sum(Y)==", sum(Y), "but n==", n, "\n")
  list(Y = Y, D = D)
  }

fExpand1 <- function(n, tp, fp, fn, tn, threshold)
  # Given scalar true/false pos/neg tests at a single threshold, deduce the number who fell above/below threshold
  {
  Y <- c( rep(1, fn + tn), rep(2, fp + tp) )
  D <- c( rep(1, fn), rep(0, tn), rep(1, tp), rep(0, fp) )
  stopifnot(length(Y) == n)
  stopifnot(length(D) == n)
  Lower = c( rep(-Inf, fn + tn), rep(threshold, fp + tp) )
  Upper = c( rep(threshold, fn + tn), rep(Inf, fp + tp) )
  list(Y = Y, D = D, Lower = Lower, Upper = Upper, n = n)
  }

fExpanddata <- function(n, vTp, vFp, vFn, vTn, vThreshold)
  # Given vectors of true/false pos/neg tests at each threshold, expand the data into patient-level data
  # Generalisation of fExpand1() for the case of more than one threshold.
  {
  # Main data
  counts <- fCategorise(n, vTp, vFp, vFn, vTn)
  Y <- rep(1, counts$Y[1])
  D <- c( rep(1, counts$D[1]), rep(0, counts$Y[1] - counts$D[1]) )
  for (j in 2:length(counts$Y)){
    Y <- c(Y, rep(j, counts$Y[j]))
    D <- c(D, rep(1, counts$D[j]), rep(0, counts$Y[j] - counts$D[j]) )
    }
  stopifnot(length(Y) == length(D))
  # Relevant thresholds
  k <- length(counts$Y)
  Lower <- rep(-Inf, counts$Y[1])
  for (j in 2:k) Lower <- c(Lower, rep(vThreshold[j - 1], counts$Y[j]))
  Upper <- rep(Inf, counts$Y[k])
  for (j in (k - 1):1) Upper <- c(rep(vThreshold[j], counts$Y[j]), Upper)
  
  list(Y = Y, D = D, Lower = Lower, Upper = Upper, n = n)
  }

fContribution <- function(i, intercept, alpha, beta, D, Lower, Upper)
  # Contribution of one observation to the loglikelihood
  {
  thetajk <- Upper[i]
  Dik <- D[i]
  numerator <- (thetajk - beta * Dik)
  denom <- exp(intercept + alpha * Dik)
  FirstTerm <- fInvlogit(numerator / denom)
  
  thetajminus1k <- Lower[i]
  Dik <- D[i]
  numerator <- (thetajminus1k - beta*Dik)
  denom <- exp(intercept+alpha*Dik)
  SecondTerm <- fInvlogit(numerator/denom)
  
  log(FirstTerm - SecondTerm)
  }

fLoglikelihood <- function(parameter)
  {
  alpha <- parameter[1] # alpha location parameter
  beta <- parameter[2] # beta scale parameter
  intercept <- parameter[3] # theta thresholds
  
  n <- length(dukic$Y)
  answer <- rep(NA, n)
  for (i in 1:n){
    answer[i] <- fContribution(i, intercept, alpha, beta, dukic$D, dukic$Lower, dukic$Upper)
    }
  sum(answer)
  }

fLoglikelihoodfixedalpha <- function(beta)
  # Loglikelihood with alpha fixed at 0
  { fLoglikelihood(c(0, beta)) }

fDukic <- function(data)
  # Make the data ready for the loglikelihood and subroutine
  {
  data <- data[order(data$study, data$cut.off), ]
  data$study <- factor(data$study)				# Strip out unused levels
  data$cut.off <- as.double(as.character(data$cut.off))	# Convert from factor to numeric
  vsStudynames <- levels(data$study)			# Identify studies in this data 
  viStudy <- c()
  iStudy <- 0
  vY <- c()
  vD <- c()
  vLower <- c()
  vUpper <- c()
  for (sStudy in vsStudynames){
    iStudy <- iStudy + 1				# Give each study an id number
    thisdata <- data[which(data$study == sStudy),]	# Extract data for this study
    iRows <- length(thisdata$n)			# Number of rows in this study (== number of thresholds)
    if (iRows > 1)	{
      lExpanded <- fExpanddata( thisdata$n[1], thisdata$TP, thisdata$FP, thisdata$FN, thisdata$TN, thisdata$cut.off )
    } else {
      lExpanded <-fExpand1( thisdata$n, thisdata$TP, thisdata$FP, thisdata$FN, thisdata$TN , thisdata$cut.off)
    }
    viStudy <- c( viStudy, rep( iStudy, lExpanded$n ) )
    vY <- c( vY, lExpanded$Y )
    vD <- c( vD, lExpanded$D )
    vLower <- c( vLower, lExpanded$Lower )
    vUpper <- c( vUpper, lExpanded$Upper )
    stopifnot(length(vY) == length(vD))
    stopifnot(length(vY) == length(viStudy))
    stopifnot(length(vY) == length(vLower))
    stopifnot(length(vY) == length(vUpper))
  }
  list(Y = vY, D = vD, Lower = vLower, Upper = vUpper, study = viStudy)
}

# ===================== SROC Simulation and Plotting Functions ===================== #

fTpr <- function(theta, intercepthat, alphahat, betahat)
  { 1 - fInvlogit( (theta - betahat) / exp(intercepthat + alphahat) ) }

fFpr <- function(theta, intercepthat)
  { 1 - fInvlogit(theta / exp(intercepthat)) }

fThetatospec <- function(theta, intercept, alpha, beta)
  { fInvlogit(theta / exp(intercept)) }

fThetatosens <- function(theta, intercept, alpha, beta)
  { 1 - fInvlogit( (theta-beta)/exp(intercept+alpha) ) }

fSpectotheta <- function(spec, intercept)
  { fLogit(spec) * exp(intercept) }
stopifnot( all.equal( 0.1, fSpectotheta( fThetatospec(0.1, 0.2), 0.2) , 1e-8  ))
stopifnot( all.equal( 0.75, fThetatospec( fSpectotheta(0.75, -0.5), -0.5 ) ))

fSenstotheta <- function(sens, intercept, alpha, beta)
  { 
  tboverea <- fLogit(1 - sens)
  theta <- tboverea * exp(intercept + alpha) + beta
  theta
  }
stopifnot( all.equal( 0.1, fSenstotheta( fThetatosens(0.1, 0.2, 0.3, 0.4), 0.2, 0.3, 0.4) , 1e-8  ))
stopifnot( all.equal( 0.75, fThetatosens( fSenstotheta(0.75, -0.5, 0.2, 0.6), -0.5, 0.2, 0.6 ), 1e-8 ))
stopifnot( all.equal( 0.1, fSpectotheta( fThetatospec(0.1, 0.2), 0.2) , 1e-8  ))
stopifnot( all.equal( 0.75, fThetatospec( fSpectotheta(0.75, -0.5), -0.5 ), 1e-8 ))

fSpectosens <- function(spec, intercept, alpha, beta)
  { 
  theta <- fSpectotheta( spec, intercept )
  sens <- fTpr(theta, intercept, alpha, beta)
  sens
  }

fSenstospec <- function(sens, intercept, alpha, beta)
  {
  tboverea <- fLogit(1 - sens)
  theta <- tboverea * exp(intercept + alpha) + beta
  spec <- 1 - fFpr(theta, intercept)
  spec
  }
stopifnot( all.equal( 0.33, fSpectosens( fSenstospec( 0.33, 0.1, 0.2, 0.3 ), 0.1, 0.2, 0.3 ), 1e-8 ))

# ===================== SROC Confidence Interval Simulation and Plotting Functions ===================== #

fSimulateCI <- function(theta, parameters, hessian, n = 1000)
  # Function for simulating confidence intervals
  {
  alpha <- parameters[1]
  beta <- parameters[2]
  intercept <- parameters[3]
  var<- -solve(hessian)
  x <- mvrnorm(n, c(alpha,beta,intercept), -solve(hessian))
  vLrminus <- vLrplus <- vSpec <- vSens <- rep(NA, n)	# Empty vectors to be filled in by the following loop:
  for (i in 1:n){
    vSens[i] <- fThetatosens(theta, x[i,3], x[i,1], x[i,2])
    vSpec[i] <- fThetatospec(theta, x[i,3], x[i,1], x[i,2])
    vLrplus[i] <- vSens[i] / (1 - vSpec[i])
    vLrminus[i] <- (1 - vSens[i]) / vSpec[i]
    }
  sensmle <- fThetatosens(theta, intercept, alpha, beta)
  specmle <- fThetatospec(theta, intercept, alpha, beta)
  list(
    theta = theta,
    lSens = list( mle = sensmle, mean = mean(vSens), lower = quantile(vSens, .025), upper = quantile(vSens, .975) ),
    lSpec = list( mle = specmle, mean = mean(vSpec), lower = quantile(vSpec, .025), upper = quantile(vSpec, .975) ),
    lLrplus = list( mle = sensmle / (1 - specmle), mean = mean(vLrplus), lower = quantile(vLrplus, .025), upper = quantile(vLrplus, .975) ),
    lLrminus = list( mle = (1 - sensmle) / specmle, mean = mean(vLrminus), lower = quantile(vLrminus, .025), upper = quantile(vLrminus, .975) )
  )
  }

fPlotsimci <- function(lSim, colour)
  # Function for plotting simulated confidence intervals
  {
  lSens <- lSim$lSens
  lSpec <- lSim$lSpec
  x <- c(1 - lSpec$mle, 1-lSpec$mle)
  y <- c(lSens$lower, lSens$upper)
  lines(x, y, col = colour, lty = "dashed")
  list(x = x, y = x)
  }

# ===================== AUC Function ===================== #

fAUC <- function(fun, a, b, tol = 1e-8, method = "simpson", verbose = TRUE)
{
  # Application of Simson's Method (or Trapezoidal rule) via Adaptive Quadrature
  # numerical integral of fun from a to b, assume a < b 
  # with tolerance tol
  n <- 4
  h <- (b - a)/4
  x <- seq(a, b, by = h)
  y <- fun(x)
  yIntegrate_internal <- function(y, h, n, method) {
    if (method == "simpson") {
      s <- y[1] + y[n + 1] + 4*sum(y[seq(2, n, by = 2)]) + 2*sum(y[seq(3, n - 1, by = 2)])
      s <- s*h/3
    } else if (method == "trapezoidal") {
      s <- h * (y[1]/2 + sum(y[2:n]) + y[n + 1]/2)
    } else {
    }		
    return(s)
  }
  
  s <- yIntegrate_internal(y, h, n, method)
  s.diff <- tol + 1 # ensures to loop at once.
  while (s.diff > tol ) {
    s.old <- s
    n <- 2*n
    h <- h/2
    y[seq(1, n + 1, by = 2)] <- y ##reuse old fun values
    y[seq(2, n, by = 2)] <- sapply(seq(a + h, b - h, by = 2*h), fun)
    s <- yIntegrate_internal(y, h, n, method)
    s.diff <- abs(s - s.old)
  }
  if (verbose) {
    cat("partition size", n, "\n")
  }
  return(s)
}

# ===================== Tabulation Functions ===================== #

fTabulatesimci <- function(lSim, sOutfilename = "")
  # Function to output a table of confidence intervals
  {
  lSens <- lSim$lSens
  lSpec <- lSim$lSpec
  lLrplus <- lSim$lLrplus
  lLrminus <- lSim$lLrminus
  cat(formatC(lSim$theta, digits = 2, format = "f"), '\t', sep = "", file = sOutfilename)
  cat(formatC(lSens$mle, digits = 2, format = "f"), " (", formatC(lSens$lower, digits = 2, format = "f"), "-", formatC(lSens$upper, digits = 2, format = "f"), ")\t", sep = "", file = sOutfilename)
  cat(formatC(lSpec$mle, digits = 2, format = "f"), " (", formatC(lSpec$lower, digits = 2, format = "f"), "-", formatC(lSpec$upper, digits = 2, format = "f"), ")\t", sep="", file = sOutfilename)
  cat(formatC(lLrplus$mle, digits = 2, format = "f"), " (", formatC(lLrplus$lower, digits = 2, format = "f"), "-", formatC(lLrplus$upper, digits = 2, format = "f"), ")\t", sep = "", file = sOutfilename)
  cat(formatC(lLrminus$mle, digits = 2, format = "f"), " (", formatC(lLrminus$lower, digits = 2, format = "f"), "-", formatC(lLrminus$upper, digits = 2, format = "f"), ")\t", sep = "", file = sOutfilename)
  cat('\n')
  }

fMyformat <- function(x)
  { 
  if (x < 1){
    formatC( x, digits = 2, format = "f") 
    } else {
    formatC( x, digits = 1, format = "f")
    }
  }

fTabulate <- function(lSimresults, cData, sOutfilename = "")
  # aRelevantdata is the subset of data rows that have the correct threshold
  {
  cat(fMyformat(lSimresults$mle), " (", fMyformat(lSimresults$lower), "-", fMyformat(lSimresults$upper), ") ", sep = "", file = sOutfilename)
  cat( "[", fMyformat(min(cData)), "-", fMyformat(max(cData)), ", n=", length(cData), "]\t", sep = "", file = sOutfilename )
  cat(fMyformat(lSimresults$mle), " (", fMyformat(lSimresults$lower), "-", fMyformat(lSimresults$upper), ")\t", sep = "", file = sOutfilename)
  }