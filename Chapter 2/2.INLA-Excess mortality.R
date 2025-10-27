### The variables in this analysis were structured as follows:
### PERIOD: time intervals (1 = 1 Jan - 30 Jun 2020, 2 = 1 July - 31 Dec 2020)
### CD_REFNIS: municipality code according to StatBel
### CD_SEX: 1 = male, 2 = female
### AGE_GROUP: age categories (1 = <25 years, 2 = 25-44 years, 3 = 45-64 years, 4 = 65-84 years, 5 = ≥85 years)
### cov_ir_log: logarithm of COVID-19 incidence rate
### N_1719: average number of deaths in 2017-2019
### idarea1: translated CD_REFNIS to be passed in INLA argument (1,...,581)

### Main libraries
library(spdep)
library(INLA)

### Create adjacency map
map_Q_nb <- poly2nb(map)
nb2INLA("MAP.adj", map_Q_nb)
gMAP <- inla.read.graph(filename = "MAP.adj")

### Model fitting with default prior
fitted.bym.def <- list()
fixed.bym.def <- list()
res.bym.def <- list()
exc.bym.def <- list()
time.bym.def <- list()

for (i in 1:2){
  data0 <- mort_20_comp%>%
    filter(PERIOD == i & AGE_GROUP %in% c("3","4","5"))
  
  start.time <- Sys.time()
  res.bym.def[[i]]  <- inla(N_20 ~ as.factor(AGE_GROUP) + as.factor(CD_SEX) + cov_ir_log +
                              f(idarea1, model = "bym", graph = gMAP, scale.model = TRUE), 
                            family = "poisson", data = data0, E = N_1719,
                            control.family = list(link = "log"),
                            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                   return.marginals.predictor=TRUE),
                            control.predictor = list(compute = TRUE))
  end.time <- Sys.time()
  time.bym.def[[paste(i)]] <- end.time - start.time
  
  fitted.bym.def[[i]] <- cbind(data0,res.bym.def[[i]]$summary.fitted.values)%>%
    mutate(PERIOD = i)
  
  fixed.bym.def[[i]] <- cbind(variable = rownames(res.bym.def[[i]]$summary.fixed),
                              data.frame(res.bym.def[[i]]$summary.fixed[,c(1:5)],
                                         row.names=NULL))%>%
    mutate(PERIOD = i)
  
  exc.bym.def[[i]] <- cbind(data0,
                            exc = sapply(res.bym.def[[i]]$marginals.fitted.values,
                                         FUN = function(marg){1 - inla.pmarginal(q = 1.2, marginal = marg)}))%>%
    mutate(PERIOD = i)
}

fitted.bym.def.all <- do.call(rbind.data.frame,fitted.bym.def)
fixed.bym.def.all <- do.call(rbind.data.frame,fixed.bym.def)
exc.bym.def.all <- do.call(rbind.data.frame,exc.bym.def)
time.bym.def.all <- do.call(rbind.data.frame,time.bym.def)

### Model fitting with PC prior
fitted.bym.pc <- list()
fixed.bym.pc <- list()
res.bym.pc <- list()
exc.bym.pc <- list()
time.bym.pc <- list()

prior <- list(
  prec = list(
    prior = "pc.prec",
    param = c(1, 0.1)))

for (i in 1:2){
  data0 <- mort_20_comp%>%
    filter(PERIOD == i & AGE_GROUP %in% c("3","4","5"))
  
  start.time <- Sys.time()
  res.bym.pc[[i]]  <- inla(N_20 ~ as.factor(AGE_GROUP) + as.factor(CD_SEX) +  cov_ir_log +
                             f(idarea1, model = "bym", graph = gMAP,
                               scale.model = TRUE,
                               hyper = list(theta1 = list(prior = "pc.prec", param = c(1, 0.1)),
                                            theta2 = list(prior = "pc.prec", param = c(1, 0.1)))), 
                           family = "poisson", data = data0, E = N_1719,
                           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                  return.marginals.predictor=TRUE),
                           control.predictor = list(compute = TRUE, link = 1))
  end.time <- Sys.time()
  time.bym.pc[[paste(i)]] <- end.time - start.time
  
  fitted.bym.pc[[i]] <- cbind(data0,res.bym.pc[[i]]$summary.fitted.values)%>%
    mutate(PERIOD = i)
  
  fixed.bym.pc[[i]] <- cbind(variable = rownames(res.bym.pc[[i]]$summary.fixed),
                             data.frame(res.bym.pc[[i]]$summary.fixed[,c(1:5)],
                                        row.names=NULL))%>%
    mutate(PERIOD = i)
  
  exc.bym.pc[[i]] <- cbind(data0,
                           exc = sapply(res.bym.pc[[i]]$marginals.fitted.values,
                                        FUN = function(marg){1 - inla.pmarginal(q = 1.2, marginal = marg)}))%>%
    mutate(PERIOD = i)
}

fitted.bym.pc.all <- do.call(rbind.data.frame,fitted.bym.pc)
fixed.bym.pc.all <- do.call(rbind.data.frame,fixed.bym.pc)
exc.bym.pc.all <- do.call(rbind.data.frame,exc.bym.pc)
time.bym.pc.all <- do.call(rbind.data.frame,time.bym.pc)
