############### Ricker Estimation Lab I - Simulation Model ###############

# Functions are provided here for the first section.
# Load these functions, then calculations will be
# given more explicitly in further sections.
Ricker <- function(S, lnalpha, beta) {
  S*exp(lnalpha - beta*S)
}

fitRicker <- function(S, R) {
  lmy <- log(R/S)
  lmfit <- lm(lmy~S)
  lnalpha_fit <- unname(lmfit$coefficients[1])
  lnalpha_p_fit <- lnalpha_fit + (sigma(lmfit)^2)/2
  beta_fit <- unname(-lmfit$coefficients[2])
  resids <- lmfit$residuals
  fits <- lmfit$fitted.values
  return(list(lnalpha_fit=lnalpha_fit, lnalpha_p_fit=lnalpha_p_fit, beta_fit=beta_fit, resids=resids, fits=fits, sigma=sigma(lmfit)))
}

simulateSR <- function(lnalpha, beta, sigS, F1, F2, sigW, N, phi, Sgoal, sigF) {
  lnalpha_p <- lnalpha + 0.5*sigW*sigW
  Seq <- lnalpha_p/beta
  # F1 <- -log(1-hrange[1])
  # F2 <- -log(1-hrange[2])
  hrange <- 1-exp(-c(F1,F2))

  # ----- initial values ----- #
  # initial value for S: Seq minus some harvest
  S <- 900 #Seq*runif(1, 1-hrange[2], 1-hrange[1])

  # initial value for observed S
  Shat <- S*rlnorm(1, sdlog=sigS)

  # initializing all other values
  redresid <- 0
  E1R <- E2R <- R <- whiteresid <- epsF <- Rgoal <- F1t <- H <- Rhat <- lnRhatShat <- fittedR <- NA

  # recursive portion...
  for(i in 2:(N+1)) {
    E1R[i] <- S[i-1]*exp(lnalpha - beta*S[i-1])
    E2R[i] <- E1R[i]*exp(phi*redresid[i-1])
    R[i] <- E2R[i]*rlnorm(1,0,sigW)
    redresid[i] <- log(R[i]/E1R[i])
    whiteresid[i] <- log(R[i]/E2R[i])
    epsF[i] <- rnorm(1,0,sigF)
    F1t[i] <- F1*exp(epsF[i])
    Rgoal[i] <- Sgoal/exp(-F1t[i])
    S[i] <- ifelse(R[i]<Rgoal[i], R[i]*exp(-F1t[i]), Sgoal+(R[i]-Rgoal[i])*exp(-F2*exp(epsF[i])))
    Shat[i] <- S[i]*rlnorm(1, sdlog=sigS)
    H[i] <- R[i]-S[i]
    Rhat[i] <- Shat[i]+H[i]
    lnRhatShat[i] <- log(Rhat[i]/Shat[i])
  }

  return(list(S=Shat[1:N],
              R=Rhat[2:(N+1)],
              Strue=S[1:N],
              Rtrue=R[2:(N+1)]))
}

calc_Smax <- function(beta) 1/beta
calc_Seq <- function(lnalpha, beta) lnalpha/beta
calc_Smsy <- function(lnalpha, beta) calc_Seq(lnalpha, beta)*(0.5-0.07*lnalpha)
calc_Umsy <- function(lnalpha) lnalpha*(0.5-0.07*lnalpha)

summary_PlotTable <- function(simdata, lnalpha, beta, sigma, Sgoal, plot=T) {
  fits <- fitRicker(S=simdata$S, R=simdata$R)
  estimates <- c(fits$lnalpha_fit,
                 fits$beta_fit,
                 fits$sigma,
                 calc_Smax(fits$beta_fit),
                 calc_Seq(fits$lnalpha_p_fit, fits$beta_fit),
                 calc_Smsy(fits$lnalpha_p_fit, fits$beta_fit),
                 calc_Umsy(fits$lnalpha_p_fit))
  truevals <- c(lnalpha,
                beta,
                sigma,
                calc_Smax(beta),
                calc_Seq(lnalpha + sigma^2/2, beta),
                calc_Smsy(lnalpha + sigma^2/2, beta),
                calc_Umsy(lnalpha + sigma^2/2))
  error <- (estimates-truevals)/truevals
  tbl <- data.frame(Estimates=estimates, Error=error)
  tbl1 <- data.frame(Estimates=estimates, Error=paste0(100*round(error,2),"%"))
  rownames(tbl) <- c("lnalpha_hat",
                     "beta_hat",
                     "sigW_hat",
                     "Smax_hat",
                     "Seq_hat",
                     "Smsy_hat",
                     "Umsy_hat")
  rownames(tbl1) <- rownames(tbl)

  # Plotting
  if(plot) {
    parmfrow <- par("mfrow") # storing graphical parameters before changing them
    par(mfrow=c(1,2))
    plotlimits <- c(0, 1.1 * max(simdata$S, simdata$R))

    # histogram of S
    hist(simdata$S, main="Observed Number of Spawners", xlim=plotlimits, breaks=seq(0,plotlimits[2],length.out=20), col=adjustcolor(4, alpha.f=.5), xlab="")

    # histogram of R
    hist(simdata$R, main="Observed Number of Recruits", xlim=plotlimits, breaks=seq(0,plotlimits[2],length.out=20), col=adjustcolor(2, alpha.f=.5), xlab="")

    # SRR plot
    plot(simdata$Strue, simdata$Rtrue, xlim=plotlimits, ylim=plotlimits, xlab="Spawners (S)", ylab="Recruits (R)", main="SRR")
    abline(0,1)
    ests <- fitRicker(S=simdata$S, R=simdata$R)
    points(simdata$Strue, Ricker(S=simdata$Strue, lnalpha=lnalpha, beta=beta), col=4, pch=15)
    points(simdata$Strue, Ricker(S=simdata$Strue, lnalpha=ests$lnalpha_fit, beta=ests$beta_fit), col=2, pch=18)
    legend("topright",legend=c("Simulated","True","Estimated"), pch=c(1,15,18), col=c(1,4,2))

    # RSR plot
    plot(simdata$Rtrue[-length(simdata$Rtrue)], simdata$Strue[-1], xlim=plotlimits, ylim=plotlimits, ylab="Spawners (S)", xlab="Recruits (R)", main="RSR")
    abline(0,1)
    abline(h=Sgoal, lty=2, lwd=2, col=2)
    legend("topleft", legend=c("Simulated data","R=S","Escapement goal"), pch=c(1,NA,NA), lty=c(NA,1,2), lwd=c(NA,1,2), col=c(1,1,2))

    par(mfrow=parmfrow) # resetting graphical parameters
  }

  if(!plot) return(tbl)
  if(plot) return(tbl1)
}

### Parts 1-6 ###

# Entering true parameter values
lnalpha <- 1.5
beta <- 0.001
sigW <- 0.4
phi <- 0
sigS <- 0.1
F1 <- 0.1 # 0.4
Sgoal <- 500
F2 <- 0.4 # 1.6
sigF <- 0.4
N <- 100


# Simulating a dataset
simdata <- simulateSR(lnalpha=lnalpha, beta=beta, sigW=sigW, phi=phi, sigS=sigS, F1=F1, F2=F2, Sgoal=Sgoal, sigF=sigF, N=N)

# Producing a summary plot, and table of parameter estimates & relative error
summary_PlotTable(simdata, lnalpha, beta, sigW, Sgoal)


# R extra: relative error of 1000 iterations of simulation & estimation!
nreps <- 1000
errormat <- matrix(nrow=nreps, ncol=7)
for(i in 1:nreps) { # this may take a few seconds to run
  simdata1 <- simulateSR(lnalpha=lnalpha, beta=beta, sigW=sigW, phi=phi, sigS=sigS, F1=F1, F2=F2, Sgoal=Sgoal, sigF=sigF, N=N)
  summarytbl <- summary_PlotTable(simdata1, lnalpha, beta, sigW, Sgoal, plot=F)
  errormat[i,] <- summarytbl[,2]
}
colnames(errormat) <- rownames(summarytbl)

# Plotting the relative error of all estimates, for all iterations of simulation & estimation.
# Parameters and reference points that are estimated more precisely will have comparatively
# narrower distributions in the boxplot.
boxplot(errormat, ylim=c(-1.5,1.5))


### Part 7 ###

# Resetting parameters to default values...
lnalpha <- 1.5
beta <- 0.001
sigW <- 0.4
phi <- 0
sigS <- 0.1
F1 <- 0.1 # 0.4
Sgoal <- 500
F2 <- 0.4 # 1.6
sigF <- 0.4
N <- 100


############### Ricker Estimation Lab II - Point Estimates ###############

### Parts 1-3 ###
# Run the following lines to simulate a fresh dataset.
# Note: data can also be read into R from a spreadsheet.
simdata <- simulateSR(lnalpha=lnalpha, beta=beta, sigW=sigW, phi=phi, sigS=sigS, F1=F1, F2=F2, Sgoal=Sgoal, sigF=sigF, N=N)
S <- tail(simdata$S, 20)  # this takes the last 20 years of our 100-year dataset
R <- tail(simdata$R, 20)


# plotting will happen at the end of Section II

log_RS <- log(R/S)  # this calculates a new vector all at once



### Part 4 ###
# Linear regression is done in R using the lm() function, in the form lm(y~x).
# Storing the results from lm() in lm_fit creates an object that we can extract information from.
lm_fit <- lm(log_RS~S)
summary(lm_fit)  # inspect the results

# plotting will happen at the end of Section II

# Extracting the coefficients.
# Note: unname() isn't needed, but some unnecessary information is carried over from lm_fit otherwise.
lnalpha_hat <- unname(lm_fit$coefficients[1])
beta_hat <- unname(-lm_fit$coefficients[2])
# beta_hat <- NA

# Bias-corrected lnalpha_p
sigma_hat <- sigma(lm_fit)
lnalpha_p_hat <- lnalpha_hat + (sigma_hat^2)/2


### Part 5 ###
Smax_hat <- 1/beta_hat
Seq_hat <- lnalpha_p_hat/beta_hat
Smsy_hat <- Seq_hat*(0.5-0.07*lnalpha_p_hat)
Umsy_hat <- lnalpha_p_hat*(0.5-0.07*lnalpha_p_hat)
MSY_hat <- Smsy_hat*exp(lnalpha_p_hat-beta_hat*Smsy_hat)-Smsy_hat
# Smax_hat <- NA
# Seq_hat <- NA
# Smsy_hat <- NA
# Umsy_hat <- NA
# MSY_hat <- NA


### Part 6-7 ###
# Fitted values and residuals can be extracted from lm_fit.
fits <- lm_fit$fitted.values
resids <- lm_fit$residuals
# # Fitted values and residuals can be extracted from lm_fit.
# # to find them try: str(lm_fit)
# fits <- NA
# resids <- NA

# The Durbin-Watson test is available in the car package.
# If the line below doesn't work, run install.packages("car") and try again.
library(car)
durbinWatsonTest(lm_fit)


### Part 8 ###
# This uses our Ricker function from above, though it could have been calculated by hand.
# Note: R does vector calculation automatically, so it returns the Rhat vector all at once.
Rhat <- Ricker(S, lnalpha_hat, beta_hat)
# # Calculate estimated recruitment by hand.
# # Note: R does vector calculation automatically, so it returns the Rhat vector all at once.
# Rhat <- NA


### Part 9 & results ###
# Printing point estimates and compare to simulated values
lnalpha_hat; lnalpha;
beta_hat; beta;
sigma_hat; sigW;
lnalpha_p_hat; (lnalpha_p <- lnalpha + sigW^2/2);
Smax_hat; 1/beta;
Seq_hat; (S_eq <- lnalpha_p / beta);
Smsy_hat; (S_msy <- S_eq * (0.5 - 0.07 * lnalpha_p));
Umsy_hat; lnalpha_p * (0.5 - 0.07 * lnalpha_p);
MSY_hat; S_msy * exp(lnalpha_p - beta * S_msy) - S_msy;
# Printing point estimates and compare to simulated values
# Blank, students complete



# Plotting...
par(mfrow=c(2,2))  # plots will now be on a 2x2 matrix
plot(S, R, xlim=c(0,max(S,R)), ylim=c(0,max(S,R)))
abline(0, 1, lty=3)  # replacement line - arguments draw a line with y-int=0, slope=1, and dotted
curve(Ricker(x, lnalpha_hat, beta_hat), add=T)  # adding a Ricker curve using our Ricker function from above
plot(S, log_RS)
abline(lm_fit)  # regression line from lm_fit
abline(h=0, lty=3)  # horizontal line at y=0
plot(S, ylim=c(0, max(S,R)), type='l', col="red")
lines(R, col="blue")
legend("topright", legend=c("S","R"), col=c("red","blue"), lty=1)
plot(resids, type='l', main="Residuals")
abline(h=0, lty=3)  # horizontal line at y=0
par(mfrow=c(1,1))





############### Ricker Estimation Lab 2 - Quantifying Uncertainty ###############


### Parts 1-5 ###
# Other methods for bootstrapping exist in R, but the sample() function works well here.
# compare the original residuals to a single resample below
data.frame(original = resids, resampled = sample(resids, replace=T)) #resampled residuals

# The structure below is called a "for loop".
# The value of i is advanced by one (1, 2, ..., B) and R performs all the calculations
# within the { } braces for each possible value of i.
B <- 10000  # how many bootstrap replicates to do.  This is upped from 1000 because computing power is cheap!
lnalpha_boot <- lnalpha_p_boot <- beta_boot <- NA  # initializing vectors to fill in within the loop
for(i in 1:B) {
  y_boot <- fits + sample(resids, replace=T)
  lm_boot <- lm(y_boot~S)
  lnalpha_boot[i] <- unname(lm_boot$coefficients[1])
  lnalpha_p_boot[i] <- lnalpha_boot[i] + 0.5*(sigma(lm_boot))^2
  beta_boot[i] <- unname(-lm_boot$coefficients[2])
}

# Censoring the impossible!  See "When beta is Poorly Defined"...
# This creates a logical vector with TRUE for impossible values, and then removes them from the bootstrap results.
impossible <- (lnalpha_boot<0) | (beta_boot<0)  # "|" = "or"
lnalpha_boot <- lnalpha_boot[!impossible]
lnalpha_p_boot <- lnalpha_p_boot[!impossible]
beta_boot <- beta_boot[!impossible]


### Part 6 ###
# Bootstrap distributions of biological reference points
# Note: these are all automatically calculated as vectors
Smax_boot <- 1/beta_boot
Seq_boot <- lnalpha_p_boot/beta_boot
Smsy_boot <- Seq_boot*(0.5-0.07*lnalpha_p_boot)
Umsy_boot <- lnalpha_p_boot*(0.5-0.07*lnalpha_p_boot)
MSY_boot <- Smsy_boot*exp(lnalpha_p_boot-beta_boot*Smsy_boot)-Smsy_boot

# Plotting as histograms
par(mfrow=c(2,2))  # Plots will now be on a 2x2 matrix
hist(lnalpha_boot)
hist(beta_boot)
hist(Seq_boot)
hist(Smsy_boot)
par(mfrow=c(1,1))


### Part 7 ###
# The quantile function can be used to calculate 10th and 90th percentiles from the bootstap distributions.
# try it for lnalpha_boot and beta_boot
quantile(lnalpha_boot, p=c(0.1, 0.9))
quantile(beta_boot, p=c(0.1, 0.9))
# quantile(lnalpha_boot, p=c(NA, NA))
# quantile(beta_boot, p=c(NA, NA))

# ...or you can do it all at once...
bootmat <- data.frame(lnalpha_boot,
                      beta_boot,
                      Smax_boot,
                      Seq_boot,
                      Smsy_boot,
                      Umsy_boot,
                      MSY_boot)   # check out head(bootmat) to see what this is
sapply(bootmat, quantile, p=c(0.1, 0.9))
# Note: to dissect the sapply() function call:
# - bootmat: apply some function to each row element of bootmat
# - quantile: what function? the quantile() function
# - p=c(0.1, 0.9): additional argument for quantile()


### Part 8 ###
diff(quantile(lnalpha_boot, p=c(0.1, 0.9))) / median(lnalpha_boot) / 2.56
# Try this yourself for beta_boot
#diff(quantile(beta_boot, p=c(0.1, 0.9)))/median(lnalpha_boot)/2.56

# It seems useful to create a NPCV function, since it looks like we'll do it a few times.
NPCV <- function(x, conf=0.8) {
  # this could all be nested in one line of code, but it's expanded to show the calculation
  conf_bounds <- c((1-conf)/2, 1-(1-conf)/2) # CI bounds from overall confidence
  quantiles <- unname(quantile(x, p=conf_bounds)) # getting the quantiles themselves
  ci_range <- diff(quantiles)  # the same as quantiles[2] - quantiles[1]
  return(ci_range/median(x)/2.56)  # what to return!
}
NPCV(lnalpha_boot)
NPCV(lnalpha_boot, conf=0.9)
sapply(bootmat, NPCV)
sapply(bootmat, NPCV, conf=0.9)


### Part 9 ###
par(mfrow=c(1,1))
plot(S, R, xlim=c(0,max(S,R)), ylim=c(0,max(S,R)))
for(i in 1:50) {
  # adding a new curve using the Ricker() function, with each bootstrap rep of lnalpha_boot and beta_boot
  curve(Ricker(x, lnalpha_boot[i], beta_boot[i]), col=adjustcolor(2,alpha.f=.3), add=T)
}
curve(Ricker(x, lnalpha, beta), lwd=2, col=4, add=T)  # adding the overall curve again
points(S, R)





############### Ricker Estimation Lab IV - Graphical Tools for Evaluating Escapement Goals ###############


### Part 1 ###
# We will be working with matrices in this section:
# - each row corresponds to one bootstrap replicate
# - each column corresponds to one prospective escapement

S_max <- 2000     # Max value of prospective escapements
S_star <- seq(1, S_max, length.out=1000)  # Prospective escapements
# expanded as a matrix
# one row for every bootstrap and one column for every prospective escapement
S_star_mat <- matrix(S_star, nrow=length(beta_boot), ncol=length(S_star), byrow=T)  # expanded as a matrix
# check dimensions
dim(S_star_mat)

# initializing the R_star matrix, then filling it in one column at a time, using bootstrap vectors all at once
R_star <- matrix(nrow=length(beta_boot), ncol=length(S_star))
for(i in 1:length(S_star)) {
  R_star[,i] <- Ricker(S_star[i], lnalpha_p_boot, beta_boot)
}

# SY_star is expanded as a matrix below
SY_star <- R_star - matrix(S_star, nrow=length(beta_boot), ncol=length(S_star), byrow=T)

# Also expanding MSY_boot and Smsy as matrix
MSY_boot_star <- matrix(MSY_boot, nrow=length(beta_boot), ncol=length(S_star))
Smsy_boot_star <- matrix(Smsy_boot, nrow=length(beta_boot), ncol=length(S_star))

# This is analogous to (c)-(d) and returns a vector of averages.
# colMeans() and rowMeans() return the means of each column or row of a matrix
# Note: we're taking the averages of ones and zeroes, just like the spreadsheet
OYP <- colMeans(SY_star >= 0.9*MSY_boot_star)  # Optimal Yield Profile
OFP <- colMeans((SY_star < 0.9*MSY_boot_star) & (S_star_mat < Smsy_boot_star))  # Overfishing Profile

# Starting a plot...
make_OYP <- function() {  # this is a shortcut for creating the plot
  plot(S_star, OYP, type='l', col="green", ylim=0:1, ylab="Probability")
  lines(S_star, OFP, col="red")
  grid()
}
par(mfrow=c(2,1))
make_OYP()


### Part 2 ###
# extracting quantiles from expected yield
quants <- c(0.05, 0.1, 0.5, 0.9, 0.95)  # which quantiles to extract
# This step pulls all quantiles at once.  To dissect the apply() function call:
# - SY_star: the SY_star matrix is what we want to apply a function to
# - 2: we want one result for each column (1=rows, 2=columns) ...this can handle higher-dimension arrays if needed
# - quantile: quantile() is the function to apply
# - p=quants: an additional argument to the quantile() function
SY_quantiles <- apply(SY_star, 2, quantile, p=quants)

# This is the first 5 columnns it returned.  There's a set of yield quantiles associated with each value of S_star.
SY_quantiles[,1:5]

make_EYP <- function() { # this is a shortcut for creating the plot
  # Making a blank plot to add lines to
  plot(NA, xlab="S", ylab="Expected Yield", xlim=range(S_star), ylim=c(0,max(SY_quantiles)))
  ltys <- c(3,2,1,2,3)  # the line type for plotting each line
  for(i in 1:5) {
    lines(S_star, SY_quantiles[i,], lty=ltys[i])
  }
  grid()
  legend("topright", legend=c("median", "80% intvl", "90% intvl"), lty=1:3)
}
make_EYP()


### Part 3 ###
# look at the OYP and EYP to come up with an escapement goal range
EG <- c(NA, NA)

par(mfrow=c(2,1))
make_OYP()
abline(v=EG, lwd=3, col="red")
make_EYP()
abline(v=EG, lwd=3, col="red")

# probabilities associated with EG endpoints
OYP[S_star %in% floor(EG)]

# yield associated with EG range
round(c(min(SY_quantiles[1, S_star>=EG[1] & S_star<=EG[2]]),
        max(SY_quantiles[5, S_star>=EG[1] & S_star<=EG[2]])))




############### Ricker Estimation Lab V - Percentile Goals ###############
# Calculate the contrast in S
max(S) / min(S)
percentiles <- c(NA, NA)  # enter your chosen percentiles here!
percEG <- round(quantile(S, percentiles))

# return the resulting EG
percEG

# and plot
par(mfrow=c(2,1))
make_OYP()
abline(v=percEG, lwd=3, col="red")
make_EYP()
abline(v=percEG, lwd=3, col="red")
par(mfrow=c(1,1))

# probabilities associated with EG endpoints
OYP[S_star %in% floor(percEG)]

# yield associated with EG range
round(c(min(SY_quantiles[1, S_star>=percEG[1] & S_star<=percEG[2]]),
        max(SY_quantiles[5, S_star>=percEG[1] & S_star<=percEG[2]])))

