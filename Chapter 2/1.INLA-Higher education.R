### The variables in this analysis were structured as follows:
### PERIOD: time intervals (1 = NA, 2 = 1-15 Jun, 3 = 16-30 Jun, 4 = 1-15 Jul,
###                         5 = 16-31 Jul, 6 = 1-15 Aug, 7 = 16-31 Aug, 8 = 1-15 Sep,
###                         9 = "16-30 Sep, 10 = 1-15 Oct, 11 = 16-31 Oct, 12 = 1-15 Nov,
###                         13 = 16-30 Nov, 14 = 1-15 Dec, 15 = 16-31 Dec)
### CD_REFNIS: municipality code according to StatBel
### C_REGIO: regional code according to StatBel (02000 = Flanders, 03000 = Walloon, 04000 = Brussels)
### AGE_CODE: age categories (≤17 years, 18-29 years, 30-39 years, 40-49 years,
###                           50-59 years, 60-69 years, 70-79 years, ≥80 years)
### SEX_CODE: M = male, F = female
### stu_ratio_centr: centered proportion of students
### gem_ink_centr: centered median income
### pop_dens_log: logarithm of population density
### pop: number of inhabitants
### idarea: translated CD_REFNIS to be passed in INLA argument (1,...,581)

### Main libraries
library(spdep)
library(INLA)

### Create adjacency matrix
map_VL <- map_BE%>%
  filter(C_REGIO == "02000")%>%
  arrange(CD_REFNIS)%>%
  mutate(idarea = 1:n())
nbVL <- poly2nb(map_VL)
nb2INLA("VL.adj", nbVL)
gVL <- inla.read.graph(filename = "VL.adj")

map_WLBRU <- map_BE%>%
  filter(C_REGIO != "02000")%>%
  arrange(CD_REFNIS)%>%
  mutate(idarea = 1:n())
nbWLBRU <- poly2nb(map_WLBRU)
nb2INLA("WLBRU.adj", nbWLBRU)
gWLBRU <- inla.read.graph(filename = "WLBRU.adj")

### Model and results
res.VL <- list()
fixed.VL <- list()
fitted.VL <- list()
ranef.VL <- list()
time.VL <- list()

for(i in c(2:15)){
  data0 <- left_join(full_2020%>%filter(C_REGIO == "02000" & PERIOD == i), 
                     map_VL%>%dplyr::select(CD_REFNIS, idarea),
                     by = "CD_REFNIS")
  start.time <- Sys.time()
  res.VL[[paste(i)]] <- inla(N ~ -1 + AGE_CODE + stu_ratio_centr:AGE_CODE + 
                               SEX_CODE + gem_ink_centr + pop_dens_log +
                               f(idarea, model = "bym", graph = gVL),
                             family = "poisson", data = data0, 
                             offset = log(pop), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                             control.predictor = list(compute = TRUE, link = 1))
  end.time <- Sys.time()
  time.VL[[paste(i)]] <- end.time - start.time
  
  fixed.VL[[paste(i)]] <- bind_cols(variable = rownames(res.VL[[paste(i)]]$summary.fixed),
                                    data.frame(res.VL[[paste(i)]]$summary.fixed[,c(1:5)],
                                               row.names = NULL))%>%
    bind_rows(., bind_cols(variable = rownames(res.VL[[paste(i)]]$summary.hyperpar)[2],
                           data.frame(res.VL[[paste(i)]]$summary.hyperpar[2,c(1:5)],
                                      row.names = NULL)))%>%
    bind_rows(., bind_cols(variable = rownames(res.VL[[paste(i)]]$summary.hyperpar)[1],
                           data.frame(res.VL[[paste(i)]]$summary.hyperpar[1,c(1:5)],
                                      row.names = NULL)))%>%
    mutate(PERIOD = i, REGIO = "Flanders")
  
  fitted.VL[[paste(i)]] <- bind_cols(data0, res.VL[[paste(i)]]$summary.fitted.values)%>%
    mutate(PERIOD = i)
  
  ranef.VL[[paste(i)]] <- left_join(map_VL, 
                                    res.VL[[paste(i)]]$summary.random$idarea%>%
                                      filter(row_number() %in% c(1:nrow(map_VL)))%>%
                                      dplyr::select(ID, mean)%>%
                                      rename(re.ss = mean),
                                    by = c("idarea" = "ID"))%>%
    mutate(idarea2 = idarea + nrow(map_VL))%>%
    left_join(., res.VL[[paste(i)]]$summary.random$idarea%>%
                filter(row_number() %in% c(nrow(map_VL)+(1:nrow(map_VL))))%>%
                dplyr::select(ID, mean)%>%
                rename(re.su = mean),
              by = c("idarea2" = "ID"))%>%
    mutate(re.total = re.ss + re.su,
           PERIOD = i)
}
fitted.VL.all <- do.call(rbind.data.frame, fitted.VL)
fixed.VL.all <- do.call(rbind.data.frame,fixed.VL)
ranef.VL.all <- do.call(rbind.data.frame,ranef.VL)
time.VL.all <- do.call(rbind.data.frame,time.VL)

res.WLBRU <- list()
fixed.WLBRU <- list()
fitted.WLBRU <- list()
ranef.WLBRU <- list()
time.WLBRU <- list()

for(i in c(2:15)){
  data0 <- left_join(full_2020%>%filter(C_REGIO != "02000" & PERIOD == i), 
                     map_WLBRU%>%dplyr::select(CD_REFNIS, idarea),
                     by = "CD_REFNIS")
  
  start.time <- Sys.time()
  res.WLBRU[[paste(i)]] <- inla(N ~ -1 + AGE_CODE + stu_ratio_centr:AGE_CODE + 
                                  SEX_CODE + gem_ink_centr + pop_dens_log +
                                  f(idarea, model = "bym", graph = gWLBRU),
                                family = "poisson", data = data0, 
                                offset = log(pop), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                                control.predictor = list(compute = TRUE, link = 1))
  end.time <- Sys.time()
  time.WLBRU[[paste(i)]] <- end.time - start.time
  
  fixed.WLBRU[[paste(i)]] <- bind_cols(variable = rownames(res.WLBRU[[paste(i)]]$summary.fixed),
                                       data.frame(res.WLBRU[[paste(i)]]$summary.fixed[,c(1:5)],
                                                  row.names = NULL))%>%
    bind_rows(., bind_cols(variable = rownames(res.WLBRU[[paste(i)]]$summary.hyperpar)[2],
                           data.frame(res.WLBRU[[paste(i)]]$summary.hyperpar[2,c(1:5)],
                                      row.names = NULL)))%>%
    bind_rows(., bind_cols(variable = rownames(res.WLBRU[[paste(i)]]$summary.hyperpar)[1],
                           data.frame(res.WLBRU[[paste(i)]]$summary.hyperpar[1,c(1:5)],
                                      row.names = NULL)))%>%
    mutate(PERIOD = i, REGIO = "Wallo-Brux")
  
  fitted.WLBRU[[paste(i)]] <- bind_cols(data0, res.WLBRU[[paste(i)]]$summary.fitted.values)%>%
    mutate(PERIOD = i)
  
  ranef.WLBRU[[paste(i)]] <- left_join(map_WLBRU, 
                                       res.WLBRU[[paste(i)]]$summary.random$idarea%>%
                                         filter(row_number() %in% c(1:nrow(map_WLBRU)))%>%
                                         dplyr::select(ID, mean)%>%
                                         rename(re.ss = mean),
                                       by = c("idarea" = "ID"))%>%
    mutate(idarea2 = idarea + nrow(map_WLBRU))%>%
    left_join(., res.WLBRU[[paste(i)]]$summary.random$idarea%>%
                filter(row_number() %in% c(nrow(map_WLBRU)+(1:nrow(map_WLBRU))))%>%
                dplyr::select(ID, mean)%>%
                rename(re.su = mean),
              by = c("idarea2" = "ID"))%>%
    mutate(re.total = re.ss + re.su,
           PERIOD = i)
}
fitted.WLBRU.all <- do.call(rbind.data.frame, fitted.WLBRU)
fixed.WLBRU.all <- do.call(rbind.data.frame,fixed.WLBRU)
ranef.WLBRU.all <- do.call(rbind.data.frame,ranef.WLBRU)
time.WLBRU.all <- do.call(rbind.data.frame,time.WLBRU)