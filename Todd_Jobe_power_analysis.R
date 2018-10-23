Model <- function(x = cobblebars,
                  type = c("normal","log","logit")){
  ## Transforms
  if (type[1] == "log")
    x$prop.woody <- log(x$prop.woody)
  else if (type[1] == "logit")
    x$prop.woody <- log(x$prop.woody / (1 - x$prop.woody))
  
  mod <- lme(prop.woody ~ year,
             data = x,
             random = ~ 1 | cobblebar/transect,
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

hps=GetHyperparam(allC20LU.lme)

hps
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

testpow=dt.power(eucC20.lme)

factoredDesign <- function(Elevs = 0.2/c(1,5,10,20),
                           Nlevs = seq(2, 10, by = 2),
                           Jlevs = seq(4, 10, by = 2),
                           Klevs = c(3, 5, 7), ...){
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
  
  return(data.frame(E, N, J, K))
}

fake <- function(N = 2, J = 6, K = 2, b = NULL, m.orig = mod,
                 transform = TRUE, ...){
  ## Simulated Data for power analysis
  ## N = Number of years
  ## J = Number of sites
  ## K = Number of stands within sites
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
      y[y > 1] <- 1 # Fix any boundary problems.
      y[y < 0] <- 0
    }
    else if (m.orig$type == "log"){
      y <- exp(eta)
      y[y > 1] <- 1
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