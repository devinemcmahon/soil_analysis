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
mod=Model(euc2deps[euc2deps$element=='C',],type='normal')

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

fixef(allC20LU.lme)
VarCorr(allC20LU.lme)
hps=GetHyperparam(allC20LU.lme)

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

dt.power(eucN,b=.25) #.941--easily detect largest changes
dt.power(eucN,b=.025) # .061
dt.power(eucN) #.16
# 15-20% chance of detecting a real change on the order of that expected
# also 5% chance of detecting a fake change

summary(budgets$Budget[budgets$Nutrient=='Ca'])
eucCa=Model(euc2deps[euc2deps$element=='Ca2',],type='log')

dt.power(eucCa,b=.4) #(median change)
dt.power(eucCa,b=.04)


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
