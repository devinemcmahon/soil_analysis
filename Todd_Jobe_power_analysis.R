# Source: 
#https://toddjobe.blogspot.com/2009/09/power-analysis-for-mixed-effect-models.html

Model <- function(x = cobblebars,
                  type = c("normal","log","logit")){
  ## Transforms
  if (type[1] == "log")
    x$stock20 <- log(x$stock20)
  else if (type[1] == "logit")
    x$stock20 <- log(x$stock20 / (1 - x$stock20))
  
  mod <- lme(stock20 ~ year,
             data = x,
             random = ~ 1 | site/stand,
             na.action = na.omit,
             control = lmeControl(opt = "optim",
                                  maxIter = 800, msMaxIter = 800)
  )
  mod$type <- type[1]
  
  return(mod)
}

# test a simple model:
eucC20.lme=lme(stock20~year,random=~1|site/stand,
                 data=euc2deps[euc2deps$element=='C',],na.action = na.omit)
#mod=Model(euc2deps[euc2deps$element=='C',],type='normal')

GetHyperparam<-function(x,b=NULL){
  ## Get the hyperparameters from the mixed effect model
  fe <- fixef(x)
  
  if(is.null(b))
    b<-fe[2] # use the data effect size if not supplied
  
  mu.a <- fe[1] 
  
  vc <- VarCorr(x)
  sigma.y <- as.numeric(vc[5, 2]) # Residual StdDev
  sigma.a <- as.numeric(vc[2, 2]) # site StdDev
  sigma.g <- as.numeric(vc[4, 2]) # site:stand StdDev
  
  hp<-c(b, mu.a, sigma.y, sigma.a, sigma.g)
  names(hp)<-c("b", "mu.a", "sigma.y", "sigma.a", "sigma.g")
  return(hp)
}

fixef(eucC20.lme)
VarCorr(eucC20.lme)
hps=GetHyperparam(eucC20.lme)

fakeModWithRestarts <- function(m.o, n = 100,  ...){
  ## A Fake Model
  withCallingHandlers({
    i <- 0
    mod <- NULL
    while (i < n & is.null(mod)){
      mod <- withRestarts({
        f <- fake(m.orig = m.o, transform = F, ...)
        return(update(m.o, data = f))
      },
      rs = function(){
        i <<- i + 1
        return(NULL)
      })
    }
    if(is.null(mod))
      warning("ExceededIterations")
    return(mod)
  },
  error = function(e){
    invokeRestart("rs")
  },
  warning = function(w){
    if(w$message == "ExceededIterations")
      cat("\n", w$message, "\n")
    else
      invokeRestart("rs")
  })
}

fake <- function(N = 2, J = 6, K = 2, b = NULL, m.orig = mod,
                 transform = TRUE, ...){
  ## Simulated Data for power analysis
  ## N = Number of years
  ## J = Number of sites
  ## K = Number of stands within sites
  # How do I make this reflect samples within stands?
  year <- rep(0:(N-1), each = J*K)
  site <- factor(rep(rep(1:J, each = K), times = N))
  stand <- factor(rep(1:K, times = N*J))
  
  ## Simulated parameters
  hp<-GetHyperparam(x=m.orig)
  if(is.null(b))
    b <- hp['b']
  g <- rnorm(J*K, 0, hp['sigma.g'])
  a <- rnorm(J*K, hp['mu.a'] + g, hp['sigma.a'])
  
  ## Simulated responses
  eta <- rnorm(J*K*N, a + b * year, hp['sigma.y'])
  if (transform){
    if (m.orig$type == "normal"){
      y <- eta
      #y[y > 1] <- 1 # Fix any boundary problems.
      #y[y < 0] <- 0
    }
    else if (m.orig$type == "log"){
      y <- exp(eta)
      #y[y > 1] <- 1
    }
    else if (m.orig$type == "logit")
      y <- exp(eta) / (1 + exp(eta))
  }
  else{
    y <- eta
  }
  
  return(data.frame(stock20 = y, year, stand, site))
}

dt.power <- function (m, n.sims = 1000, alpha=0.05, ...){
  ## Calculate power for a particular sampling design
  signif<-rep(NA, n.sims)
  for(i in 1:n.sims){
    lme.power <- fakeModWithRestarts(m.o = m, ...)
    if(!is.null(lme.power))
      signif[i] <- summary(lme.power)$tTable[2, 5] < alpha
  }
  power <- mean(signif, na.rm = T)
  return(power)
}

testpow=dt.power(euc20.lme) #.184. so not great. (again: .187)

eucN=Model(euc2deps[euc2deps$element=='N',],type='log')
summary(budgets$Budget[budgets$Nutrient=='N'])
# N budgets (estimated) range from about -.25 to +.25 Mg ha-1
# mean is -.04
# mean change in ln(value) = ln(N16/N04) = .093
# expected change in ln(value)
tapply(stkchgs$efs20,stkchgs$element,mean)
tapply(stkchgs$budget,stkchgs$element,mean) # budget effect size
tapply(stkchgs$efs20,stkchgs$element,median) # log budget effect size
tapply(stkchgs$chgln20,stkchgs$element,median) # log observed effect size
tapply(stkchgs$chg20,stkchgs$element,median) # observed effect size

# N: .0005 (median .023)

dt.power(eucN,b=.023) # .061
GetHyperparam(eucN) #b=.07
dt.power(eucN,b=.093) #.247
dt.power(eucN,b=.0005) #0.05
dt.power(eucN) #.16
# expected effect
# 15-25% chance of detecting a real change on the order of that observed
# 5-6% change of detecting change on order of that expected
eucC=Model(euc2deps[euc2deps$element=='C',],type='normal')

median(abs(shorttstk$stock20_16[shorttstk$element=='C' & shorttstk$LU=='E']-
         shorttstk$stock20_04[shorttstk$element=='C'& shorttstk$LU=='E']))
# 6.086 w/o abs, 8.553 with
dt.power(eucC,b=6.086) # .407
dt.power(eucC) # .197. hm. Is this really testing what I want it to?
# Yes, because I'm looking for an overall effect, so opposing effects cancel
#   even if the change is large in a given stand

summary(budgets$Budget[budgets$Nutrient=='Ca'])
eucCa=Model(euc2deps[euc2deps$element=='Ca2',],type='log')

dt.power(eucCa,b=.4) #(median change)
# no, with log, median change is 1.97
dt.power(eucCa,b=1.97) #1
dt.power(eucCa,b=.04)
dt.power(eucCa)
dt.power(eucCa,b=1.883)

eucP=Model(euc2deps[euc2deps$element=='P2',])
dt.power(eucP)
dt.power(eucP,b=.039) #median observed change
dt.power(eucP,b=.075) # median budget change--should be mean? that's .070

eucK=Model(euc2deps[euc2deps$element=='K'&
                      euc2deps$site!='Bp',],type='log')
qqr(eucK) # tails still off with Bp--without Bp, nice with log
dt.power(eucK)
mean(log(tstock$stock20_16[tstock$element=='K' & tstock$LU=='E' & 
                             tstock$site!='Bp']/
           tstock$stock20_04[tstock$element=='K' & tstock$LU=='E' & 
                               tstock$site!='Bp'])) # mean .0809, med .106
dt.power(eucK,b=.092) #median observed change (ln)
dt.power(eucK,b=.232) # median budget (close to mean)
dt.power(eucK, b=.081)

eucCu=Model(euc2deps[euc2deps$element=='Cu',],type='normal')
qqr(eucCu) # some outliers; better without log?
dt.power(eucCu) # .055
summary(eucCu) # change is .0005 Mg? ha-1 i.e. half a kg
mean(tstock$stock20_16[tstock$element=='Cu' & tstock$LU=='E']-
             tstock$stock20_04[tstock$element=='Cu'& tstock$LU=='E']) 
#.0013 Mg, 13 kg = median, mean is .0006, ok
dt.power(eucCu,b=.0013) # .113 not very useful if overall effect != median

eucZn=Model(euc2deps[euc2deps$element=='Zn',],type='log')
qqr(eucZn) # tails still off
dt.power(eucZn) # .055
summary(eucZn) # 
mean(log(tstock$stock20_16[tstock$element=='Zn' & tstock$LU=='E']/
       tstock$stock20_04[tstock$element=='Zn'& tstock$LU=='E']))
#med .0607, mean.127
dt.power(eucZn,b=.061) # .054, although I didn't log-transform it


factoredDesign <- function(Elevs = 0.25/c(.5,1,2,5,10),
                           Nlevs = 2,
                           Jlevs = seq(2, 10, by = 2),
                           Klevs = c(1,2), ...){
  ## Generates factored series of sampling designs for simulation
  ## of data that follow a particular model.
  ## Inputs:
  ##   Elevs - vector of effect sizes for the slope parameter.
  ##   Nlevs - vector of number of years to sample.
  ##   Jlevs - vector of number of cobblebars to sample.
  ##   Klevs - vector of number of transects to sample.
  ## Results:
  ##   Data frame with where columns are the factors and
  ##   rows are the designs.
  
  # Level lengths
  lE <- length(Elevs)
  lN <- length(Nlevs)
  lJ <- length(Jlevs)
  lK <- length(Klevs)
  
  # Generate repeated vectors for each factor
  E <- rep(Elevs, each = lN*lJ*lK)
  N <- rep(rep(Nlevs, each = lJ*lK), times = lE)
  J <- rep(rep(Jlevs, each = lK), times = lE*lN)
  K <- rep(Klevs, times = lE*lN*lJ)
  # This isn't what I want to test--what about samples per transect?
  # His nesting scheme is not the same as mine: 
  #   more like resampling individual reps in a stand
  # Need to modify this code to estimate power with different no. samples
  # Even then, hyperparameters constrained by variances in limited sampling
  # Changing no. samples within a stand won't change the power of my lme?
  #   No, it should reduce the residual variance
  # In paper, address/focus on detection of changes within-stand (t-test)
  #   as will be useful for companies? Or maybe they want to know how multiple
  #   stands are changing on average, too.

  
  return(data.frame(E, N, J, K))
}



powerAnalysis <- function(parallel = T, ...){
  ## Full Power Analysis
  
  ## Parallel
  if(parallel){
    closeAllConnections()
    cl <- makeCluster(7, type = "SOCK")
    on.exit(closeAllConnections())
    clusterEvalQ(cl, source("cobblebars2.r"))
  }
  
  ## The simulations
  dat <- factoredDesign(...)
  
  if (parallel){
    dat$power <- parRapply(cl, dat, function(x,...){
      dt.power(N = x[2], J = x[3], K = x[4], b = x[1], ...)
    }, ...)
  } else {
    dat$power <- apply(dat, 1, function(x, ...){
      dt.power(N = x[2], J = x[3], K = x[4], b = x[1], ...)
    }, ...)
  }
  
  return(dat)
}

powsN=powerAnalysis(parallel = F,m=eucN)
# Yes, that takes ages to run
#saveRDS(powsN,'powsn.Rds')

plotPower <- function(dt){
  xyplot(power~N|J*K, data = dt, groups = E,
         panel = function(...){panel.xyplot(...)
           panel.abline(h = 0.8, lty = 2)},
         type = c("p", "l"),
         xlab = "sampling years",
         ylab = "power",
         strip = strip.custom(var.name = c("C", "T"),
                              strip.levels = c(T, T)),
         auto.key = T
  )
}
plotPower(powsN)


# other power tests
power.t.test(n=9,delta=mean(stkchgs$chg20[stkchgs$element=='N']),
             sd=mean(stkchgs$sdchg20[stkchgs$element=='N']),
             sig.level = .05,type='paired') # power=.488
power.t.test(n=9,delta=mean(stkchgs$budget[stkchgs$element=='N']),
             sd=mean(stkchgs$sdchg20[stkchgs$element=='N']),
             sig.level = .05,type='paired') # .058
