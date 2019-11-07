library(shiny)

ui <- fluidPage(

  sidebarLayout(
    sidebarPanel(
      h3("Inputs"),
      checkboxInput("showref","Show biological reference points",value=T),
      checkboxInput("sim","Simulate data"),
      actionButton("sim1","generate new points..."),

      tabsetPanel(
        tabPanel("Ricker",

                 h3("Ricker Parameters:"),
                 sliderInput("lnalpha",
                             "Productivity ln(alpha):",
                             min = 0.01,
                             max = 3,
                             value = 2),
                 sliderInput("beta",
                             "Density dependence beta:",
                             min = 0.00001,
                             max = 0.001,
                             value = 0.0005)
        ),



        tabPanel("Data Sim",
                 h3("Simulate data:"),

                 # sliderInput("N","Number of observations:", min=2, max=100, value=30, ticks=T),
                 # sliderInput("sig","Process error lognormal sigma:", min=0, max=1, value=0.3, step=.01, ticks=T),
                 # # sliderInput("Srange","Data range of S (S/Seq)", min=0, max=2, value=c(0.2,1.2), step=.01),
                 # sliderInput("hrange","Range of harvest rates:", min=0, max=1, value=c(0.5,0.9), step=.01, ticks=T),
                 # sliderInput("cvS","Measurement error cv(S):", min=0, max=0.5, value=0, step=.01, ticks=T),
                 # sliderInput("phi","Autocorrelation phi:", min=0, max=1, value=0, step=.01, ticks=T)
                 div(style="height: 70px;", sliderInput("N","Number of observations:", min=2, max=100, value=30, ticks=F)),
                 div(style="height: 70px;", sliderInput("sigW","Process error lognormal sigma:", min=0, max=1, value=0.3, step=.01, ticks=F)),
                 div(style="height: 70px;", sliderInput("sigS","Measurement error sigma_S:", min=0, max=0.5, value=0, step=.01, ticks=F)),
                 div(style="height: 70px;", sliderInput("phi","Autocorrelation phi:", min=0, max=1, value=0, step=.01, ticks=F)),
                 div(style="height: 70px;", sliderInput("hrange","Harv rates: typical and liberalized", min=0, max=1, value=c(0.25,0.75), step=.01, ticks=F)),
                 div(style="height: 70px;", sliderInput("Sgoal","S Goal for liberalization:", min=0, max=10000, value=2000, step=1, ticks=F)),
                 div(style="height: 70px;", sliderInput("sigF","Harvest lognormal sigma_F:", min=0, max=1, value=0.4, step=.01, ticks=F))
        )
      ),
      sliderInput("maxS","plot size:", min=0,max=50000,value=10000)
    ),

    mainPanel(
      h3("Outputs"),

      tabsetPanel(tabPanel("Ricker & Yield Profiles",
                           plotOutput("RickerPlot", height="400px", width="600px"),
                           plotOutput("YieldPlot", height="400px", width="600px")
      ),
      tabPanel("Fit & Residuals",
               plotOutput("ResidPlot", height="600px", width="600px")
      ),
      tabPanel("Naive Bootstrap Histograms",
               plotOutput("BootPlot", height="600px", width="600px")
      ),
      tabPanel("Meta-Simulation",
               h5("In this section, many iterations of possible datasets are simulated under given Ricker and data parameters."),
               h5("True values (blue) can be compared to the true dispersion (black, estimated as the middle 80% of point estimates),
                  and estimated dispersion (red, estimated as the average of 80% bootstrap CI endpoints)."),
               div(style="height: 70px;", sliderInput("metareps","Number of iterations:", min=10,max=1000,value=100,step=10, ticks=F)),
               div(style="height: 70px;", sliderInput("metaB","Number of bootstrap reps per iteration:", min=10,max=1000,value=100,step=10, ticks=F)),
               actionButton("runmeta","run meta-simulation..."),
               plotOutput("MetaPlot", height="500px", width="600px")
               ))
    )
  )

)

# Define server logic
server <- function(input, output) {

  # library(car)

  # ----- functions ----- #

  Ricker <- function(x, lnalpha=input$lnalpha, beta=input$beta) x*exp(lnalpha - beta*x)

  fitRicker <- function(S, R) {
    lmy <- log(R/S)
    lmfit <- lm(lmy~S)
    lnalpha_fit <- unname(lmfit$coefficients[1])
    lnalpha_p_fit <- lnalpha_fit + (sigma(lmfit)^2)/2
    beta_fit <- unname(-lmfit$coefficients[2])
    resids <- lmfit$residuals
    fits <- lmfit$fitted.values
    return(list(lnalpha_fit=lnalpha_fit, lnalpha_p_fit=lnalpha_p_fit, beta_fit=beta_fit, resids=resids, fits=fits))
  }

  bootRicker <- function(S, R, B=1000) {
    firstfit <- fitRicker(S=S, R=R)
    fits <- firstfit$fits
    resids <- firstfit$resids
    lnalpha_boot <- lnalpha_p_boot <- beta_boot <- rep(NA, B)
    for(i in 1:B) {
      lmy <- fits + sample(resids, replace=T)
      lmfit <- lm(lmy~S)
      lnalpha_boot[i] <- unname(lmfit$coefficients[1])
      lnalpha_p_boot[i] <- lnalpha_boot[i] + (sigma(lmfit)^2)/2
      beta_boot[i] <- unname(-lmfit$coefficients[2])
    }

    impossible <- lnalpha_boot<0 | beta_boot<0  # censor the impossible
    return(list(lnalpha_boot=lnalpha_boot[!impossible],
                lnalpha_p_boot=lnalpha_p_boot[!impossible],
                beta_boot=beta_boot[!impossible]))
  }

  simulateSR <- function(lnalpha, beta, sigS, hrange, sigW, N, phi, Sgoal, sigF) {             ##### add Sgoal and sigF to calls of this function
    # calculate lognormal sigma in measurement error for S from cv(S)
    # siglog <- sqrt(log(cvS^2 + 1))
    lnalpha_p <- lnalpha + 0.5*sigW*sigW
    Seq <- lnalpha_p/beta
    F1 <- -log(1-hrange[1])
    F2 <- -log(1-hrange[2])
    # hrange <- 1-exp(-c(F1,F2))

    # ----- initial values ----- #
    # initial value for S: Seq minus some harvest
    S <- Seq*runif(1, 1-hrange[2], 1-hrange[1])       ### this could be betterized, but I guess it's not so bad

    # initial value for observed S
    Shat <- S*rlnorm(1, sdlog=sigS)

    # initializing all other values
    redresid <- 0 ## should this be betterized?
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

  metaSimulate <- function(lnalpha, beta, sigS, hrange, sigW, N, phi, Sgoal, sigF, boot=T, reps=500, B=500) {
    estimates <- cilo <- cihi <- as.data.frame(matrix(nrow=reps, ncol=6))
    names(estimates) <- names(cilo) <- names(cihi) <- c("lnalpha_p","beta","Smsy","Smax","Seq","MSY")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...', value = 0, {
                   for(i in 1:reps) {
                     incProgress(i/reps)   ################
                     thesim <- simulateSR(lnalpha=lnalpha,
                                          beta=beta,
                                          sigS=sigS,
                                          hrange=hrange,
                                          sigW=sigW,
                                          N=N,
                                          phi=phi,
                                          Sgoal=Sgoal,
                                          sigF=sigF)
                     thefit <- fitRicker(S=thesim$S, R=thesim$R)
                     estimates$lnalpha_p[i] <- thefit$lnalpha_p_fit
                     estimates$beta[i] <- thefit$beta_fit
                     if(boot) {
                       theboot <- bootRicker(S=thesim$S, R=thesim$R, B=B)
                       Seq_boot <- theboot$lnalpha_p_boot/theboot$beta_boot
                       Smax_boot <- 1/theboot$beta_boot
                       Smsy_boot <- Seq_boot*(0.5-0.07*theboot$lnalpha_p_boot)
                       MSY_boot <- Ricker(Smsy_boot, theboot$lnalpha_p_boot, theboot$beta_boot) - Smsy_boot
                       cilo[i,] <- sapply(list(theboot$lnalpha_p_boot, theboot$beta_boot, Smsy_boot, Smax_boot, Seq_boot, MSY_boot),
                                          quantile, p=0.1)
                       cihi[i,] <- sapply(list(theboot$lnalpha_p_boot, theboot$beta_boot, Smsy_boot, Smax_boot, Seq_boot, MSY_boot),
                                          quantile, p=0.9)
                     }
                   }})
    estimates$Seq <- estimates$lnalpha_p/estimates$beta
    estimates$Smax <- 1/estimates$beta
    estimates$Smsy <- estimates$Seq*(0.5-0.07*estimates$lnalpha_p)
    estimates$MSY <- Ricker(estimates$Smsy, estimates$lnalpha_p, estimates$beta) - estimates$Smsy
    if(boot) out <- list(estimates=estimates, cilo=cilo, cihi=cihi)
    if(!boot) out <- estimates
    return(out)
  }

  plotMeta <- function(x, lnalpha, beta) {
    Seq <- lnalpha/beta
    Smax <- 1/beta
    Smsy <- Seq*(0.5-0.07*lnalpha)
    MSY <- Ricker(Smsy, lnalpha, beta) - Smsy
    # histt <- function(y,Y,mnci=NULL,...) {
    #     maxthing <- max(abs(range(y,mnci,na.rm=T)-Y))
    #     xlimthing <- Y+c(-1,1)*maxthing
    #     hist(y, xlim=xlimthing, xlab="", ...=...)
    #     abline(v=Y,col=4,lwd=3)
    #     abline(v=mnci, col=2, lwd=2)
    #     abline(v=quantile(y, p=c(.1,.9)))
    # }
    histt <- function(y,Y,mnci=NULL,...) {
      maxthing <- max(abs(range(quantile(y,p=c(0.1,.9)),mnci,na.rm=T)-Y))
      xlimthing <- Y+c(-1,1)*maxthing
      plot(NA, xlim=xlimthing, ylim=0:1, yaxt="n", xlab="", ...=...)
      abline(v=Y,col=4,lwd=3)
      # abline(v=mnci, col=2, lwd=2)
      lines(mnci, rep(0.3,2), col=2, lwd=2, lend=1)
      # abline(v=quantile(y, p=c(.1,.9)))
      lines(quantile(y, p=c(.1,.9)), rep(0.6,2), col=1, lwd=2, lend=1)
    }
    par(mfrow=c(6,1), mar=c(3,2,1,2))
    histt(x$estimates$lnalpha_p, lnalpha, c(mean(x$cilo$lnalpha_p), mean(x$cihi$lnalpha_p)), main="lnalpha")
    if(mean(x$estimates$lnalpha_p)>=lnalpha) legend("topleft",col=c(1,2),lwd=2,legend=c("true dispersion","est dispersion"))
    if(mean(x$estimates$lnalpha_p)<lnalpha) legend("topright",col=c(1,2),lwd=2,legend=c("true dispersion","est dispersion"))
    histt(x$estimates$beta, beta, c(mean(x$cilo$beta), mean(x$cihi$beta)), main="beta")
    histt(x$estimates$Smsy, Smsy, c(mean(x$cilo$Smsy), mean(x$cihi$Smsy)), main="Smsy")
    histt(x$estimates$Smax, Smax, c(mean(x$cilo$Smax), mean(x$cihi$Smax)), main="Smax")
    histt(x$estimates$Seq, Seq, c(mean(x$cilo$Seq), mean(x$cihi$Seq)), main="Seq")
    histt(x$estimates$MSY, MSY, c(mean(x$cilo$MSY), mean(x$cihi$MSY)), main="MSY")
  }

  plotRicker <- function(lnalpha, beta, maxS, showref=T, sim=T, S, R, lnalpha_boot, lnalpha_p_boot, beta_boot) {
    DrawRef <- function(x,...) lines(rep(x,2),c(0,Ricker(x)), ...=...)
    Seq <- lnalpha/beta
    Smax <- 1/beta
    Smsy <- Seq*(0.5-0.07*lnalpha)
    MSY <- Ricker(Smsy, lnalpha, beta) - Smsy

    plot(NA, xlim=c(0,maxS), ylim=c(-.05*maxS,maxS), xlab="S", ylab="R")
    abline(0,1,col="grey")
    alphang <- atan(exp(lnalpha))
    lines(maxS*.4*c(0,cos(alphang)), maxS*.4*c(0,sin(alphang)), lty=2, col="grey30")
    if(showref) {
      DrawRef(Smsy, col=4)
      DrawRef(Smax, col=4, lty=2)
      DrawRef(Seq, col=4, lty=3)
      text(c(Smsy,Smax,Seq),rep(-.03*input$maxS,3),labels=c("Smsy","Smax","Seq"),col=4,cex=.8) #,pos=c(2,4,4)
    }
    if(input$sim) {
      points(S, R)
      Rfit <- fitRicker(S=S, R=R)
      for(i in 1:min(length(beta_boot),100)) curve(Ricker(x, lnalpha=lnalpha_boot[i], beta=beta_boot[i]), add=T, col=adjustcolor(2,alpha.f=.25))
      curve(Ricker(x, lnalpha=Rfit$lnalpha_fit, beta=Rfit$beta_fit), add=T, col=2)
      Seq_boot <- lnalpha_p_boot/beta_boot
      Smax_boot <- 1/beta_boot
      Smsy_boot <- Seq_boot*(0.5-0.07*lnalpha_p_boot)
      if(input$showref) {
        lines(quantile(Smsy_boot, c(.025,.975)), rep(maxS*.01, 2), col=2)
        lines(quantile(Smsy_boot, c(.25,.75)), rep(maxS*.01, 2), col=2, lwd=3, lend=1)
        lines(quantile(Smax_boot, c(.025,.975)), rep(maxS*.03, 2), col=2)
        lines(quantile(Smax_boot, c(.25,.75)), rep(maxS*.03, 2), col=2, lwd=3, lend=1)
        lines(quantile(Seq_boot, c(.025,.975)), rep(maxS*.05, 2), col=2)
        lines(quantile(Seq_boot, c(.25,.75)), rep(maxS*.05, 2), col=2, lwd=3, lend=1)
      }
    }
    curve(Ricker(x), col=4, lwd=2, add=T)
  }



  # ----- reactive values ---- #

  sim <- reactive({
    trick1 <- input$sim1
    simSR <- simulateSR(lnalpha=input$lnalpha,
                        beta=input$beta,
                        sigS=input$sigS,
                        hrange=input$hrange,
                        sigW=input$sigW,
                        N=input$N,
                        phi=input$phi,
                        Sgoal=input$Sgoal,
                        sigF=input$sigF)
    return(simSR)
  })

  Rboot1 <- reactive({
    trick1 <- input$sim1
    if(input$sim) {# & input$radio>=2) {
      Rboot <- bootRicker(S=sim()$S, R=sim()$R)
    }
    Rboot
  })

  metasim <- eventReactive(input$runmeta, {
    # trick1 <- input$runmeta
    # if(input$sim) {# & input$radio>=2) {
    meta <- metaSimulate(lnalpha=input$lnalpha,
                         beta=input$beta,
                         sigS=input$sigS,
                         hrange=input$hrange,
                         sigW=input$sigW,
                         N=input$N,
                         phi=input$phi,
                         Sgoal=input$Sgoal,
                         sigF=input$sigF,
                         boot=T,
                         reps=input$metareps,
                         B=input$metaB)
    # }
    meta
  })



  # ----- reactive plots ---- #

  output$RickerPlot <- renderPlot({
    plotRicker(lnalpha=input$lnalpha,
               beta=input$beta,
               maxS=input$maxS,
               showref=input$showref,
               sim=input$sim,
               S=sim()$S,
               R=sim()$R,
               lnalpha_boot=Rboot1()$lnalpha_boot,
               lnalpha_p_boot=Rboot1()$lnalpha_p_boot,
               beta_boot=Rboot1()$beta_boot)
  })

  output$YieldPlot <- renderPlot({
    if(input$sim) {# & input$radio==3) {
      Rboot <- Rboot1()
      par(mfrow=c(3,1), mar=c(2,4,1,2)+.1)
      lyp <- 200 # length of the vector considered for yield profiles
      Syp <- seq(1, input$maxS, length.out=lyp) # S considered
      Ryp <- matrix(nrow=length(Rboot$lnalpha_boot), ncol=lyp)   # each row is a bootstrap rep, each column is an Syp
      for(i in 1:lyp) Ryp[,i] <- Ricker(Syp[i], lnalpha=Rboot$lnalpha_boot, beta=Rboot$beta_boot)
      Yyp <- Ryp-matrix(Syp, nrow=length(Rboot$lnalpha_boot), ncol=lyp, byrow=T)
      Seq_boot <- Rboot$lnalpha_p_boot/Rboot$beta_boot
      Smax_boot <- 1/Rboot$beta_boot
      Smsy_boot <- Seq_boot*(0.5-0.07*Rboot$lnalpha_p_boot)
      MSY_boot <- Ricker(Smsy_boot, lnalpha=Rboot$lnalpha_boot, beta=Rboot$beta)-Smsy_boot
      probs <- c(.9,.8,.7)

      Seq <- input$lnalpha/input$beta
      Smax <- 1/input$beta
      Smsy <- Seq*(0.5-0.07*input$lnalpha)

      addbootCIs <- function() {
        maxy <- par("usr")[4]
        thecol=adjustcolor(2, alpha.f=.4)
        lines(quantile(Smsy_boot, c(.025,.975)), rep(maxy*.04, 2), col=thecol)
        lines(quantile(Smsy_boot, c(.25,.75)), rep(maxy*.04, 2), col=thecol, lwd=3, lend=1)
        lines(quantile(Smax_boot, c(.025,.975)), rep(maxy*.08, 2), col=thecol)
        lines(quantile(Smax_boot, c(.25,.75)), rep(maxy*.08, 2), col=thecol, lwd=3, lend=1)
        lines(quantile(Seq_boot, c(.025,.975)), rep(maxy*.12, 2), col=thecol)
        lines(quantile(Seq_boot, c(.25,.75)), rep(maxy*.12, 2), col=thecol, lwd=3, lend=1)
      }

      plot(NA, xlim=c(0, input$maxS), ylim=0:1, main="Optimal Yield Profile", xaxt="n")
      grid()
      legend("topright", lty=1:length(probs), legend=paste0(100*probs,"% of MSY"), title="Probability of achieving")
      if(input$showref) {
        abline(v=c(Smsy, Smax, Seq), lty=1:3, col=adjustcolor(4, alpha.f=.4), lwd=2)
        text(Smsy, par("usr")[4], labels="true Smsy", col=4, xpd=NA, pos=3)
        addbootCIs()
      }
      for(i in 1:length(probs)) {
        pp <- colMeans(Yyp>=(MSY_boot*probs[i]))
        lines(Syp,pp,lty=i)
      }

      plot(NA, xlim=c(0, input$maxS), ylim=0:1, main="Overfishing Profile", xaxt="n")
      grid()
      legend("topright", lty=1:length(probs), legend=paste0(100*probs,"% of MSY"), title="Probability of reduction to")
      if(input$showref) {
        abline(v=c(Smsy, Smax, Seq), lty=1:3, col=adjustcolor(4, alpha.f=.4), lwd=2)
        addbootCIs()
      }
      for(i in 1:length(probs)) {
        pp <- colMeans(Yyp<=(MSY_boot*probs[i]) & outer(Smsy_boot, Syp, function(X,Y) X>=Y))
        lines(Syp,pp,lty=i)
      }

      EYyp <- apply(Yyp, 2, quantile, p=c((1-probs)/2, 1-((1-rev(probs))/2)))
      ltys <- c(length(probs):1, 1:length(probs))
      plot(NA, xlim=c(0, input$maxS), ylim=c(0,max(EYyp,na.rm=T)), main="Expected Yield Profile", xlab="S")
      grid()
      legend("topright", lty=c(length(probs):1,1), legend=c(paste0(100*probs,"% intvl"),"median"))
      if(input$showref) {
        abline(v=c(Smsy, Smax, Seq), lty=1:3, col=adjustcolor(4, alpha.f=.4), lwd=2)
        addbootCIs()
      }
      for(i in 1:(2*length(probs))) {
        lines(Syp,EYyp[i,],lty=ltys[i])
      }
      lines(Syp,apply(Yyp, 2, quantile, p=.5,lwd=3))
    }
  })

  output$ResidPlot <- renderPlot({
    if(input$sim) {
      par(mfrow=c(2,2), mar=c(4,4,2,2))
      plot(sim()$S, sim()$R, xlab="S", ylab="R", pch=16, col=adjustcolor(2,alpha.f=.4), xlim=c(0,input$maxS), ylim=c(0,input$maxS))
      if(input$sigS>0) {
        points(sim()$Strue, sim()$Rtrue, pch=16, col=adjustcolor(4,alpha.f=.4))
        segments(sim()$Strue, sim()$Rtrue, sim()$S, sim()$R, col="grey")
      }
      curve(Ricker(x, input$lnalpha, input$beta), add=T, col=4, lwd=1)
      Rfit <- fitRicker(S=sim()$S, R=sim()$R)
      # print(Rfit)
      curve(Ricker(x, lnalpha=Rfit$lnalpha_fit, beta=Rfit$beta_fit), add=T, col=2)
      legend("topright",legend=c("True","Fit"),col=c(4,2),lwd=2)

      logRS <- log(sim()$R/sim()$S)
      logRStrue <- log(sim()$Rtrue/sim()$Strue)
      plot(sim()$S, logRS, xlab="S", ylab="log(R/S)", pch=16, col=adjustcolor(2,alpha.f=.4), xlim=c(0,max(sim()$S,sim()$Strue)), ylim=range(logRS,logRStrue,na.rm=T))
      if(input$sigS>0) {
        points(sim()$Strue, logRStrue, pch=16, col=adjustcolor(4,alpha.f=.4))
        segments(sim()$Strue, logRStrue, sim()$S, logRS, col="grey")
      }
      abline(Rfit$lnalpha_fit, -Rfit$beta_fit, col=2)
      abline(input$lnalpha, -input$beta, col=4)
      legend("topright",legend=c("True","Fit"),col=c(4,2),lwd=2)

      plotmax <- max(abs(Rfit$resids))
      plot(Rfit$resids, ylim=c(-plotmax,1.2*plotmax), type='l', main="Residuals", xlab="Year", ylab="Regression residual")
      abline(h=0, lty=3)
      text(1, 1.2*plotmax, paste("D-W stat:", round(car::durbinWatsonTest(Rfit$resids),2)), pos=4)
      text(1, 1.1*plotmax, paste("Observed autocorrelation:", round(cor(Rfit$resids[-1],Rfit$resids[-input$N]),2)), pos=4)

      resids <- Rfit$resids
      acfplot <- acf(resids, plot=F)
      plot(acfplot, main="Autocorrelation Plot")
    }
  })

  output$BootPlot <- renderPlot({
    if(input$sim) {# & input$radio>=2) {

      Rboot <- Rboot1()#bootRicker(S=Ssim, R=Rsim)
      Seq_boot <- Rboot$lnalpha_p_boot/Rboot$beta_boot
      Smax_boot <- 1/Rboot$beta_boot
      Smsy_boot <- Seq_boot*(0.5-0.07*Rboot$lnalpha_p_boot)
      MSY_boot <- Ricker(Smsy_boot, Rboot$lnalpha_p_boot, Rboot$beta_boot) - Smsy_boot

      Seq <- input$lnalpha/input$beta
      Smax <- 1/input$beta
      Smsy <- Seq*(0.5-0.07*input$lnalpha)
      MSY <- Ricker(Smsy) - Smsy

      histt <- function(x,X,...) {
        maxthing <- max(abs(range(x,na.rm=T)-X))
        xlimthing <- X+c(-1,1)*maxthing
        hist(x, xlim=xlimthing, xlab="", ...=...)
        abline(v=X,col=4,lwd=3)
      }
      par(mfrow=c(6,1), mar=c(3,2,1,2))
      histt(Rboot$lnalpha_boot, input$lnalpha, main="lnalpha")
      histt(Rboot$beta_boot, input$beta, main="beta")
      histt(Smsy_boot, Smsy, main="Smsy")
      histt(Smax_boot, Smax, main="Smax")
      histt(Seq_boot, Seq, main="Seq")
      histt(MSY_boot, MSY, main="MSY")
    }
  })

  output$MetaPlot <- renderPlot({
    plotMeta(metasim(), input$lnalpha, input$beta)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
