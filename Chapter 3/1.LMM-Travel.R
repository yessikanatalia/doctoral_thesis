### The variables in this analysis were structured as follows:
### Week: weekly time interval based on ISO 8601 definition
### C_PROVI: province code according to StatBel
### IC7: weekly incidence of new cases
### IH7: weekly incidence of hospitalizations
### fullvax: vaccination coverage
### fullvaxLL: lagged vaccination coverage
### fullvaxi: inverse of vaccination coverage
### fullvaxiLL: lagged inverse of vaccination coverage
### travel: incoming travel rate
### pos_perc: positivity rate of the incoming travelers
### SI: median of the stringency index

### Main libraries
library(lme4)

### Variable definition
outcome <- c("IC7", "IH7")
inverse <- c("noinv", "inv")
lag <- as.character(c(0:24))
lag.alt <- as.character(c(0:4,10:12))

res.all <- list()
mod.all <- list()
res.20 <- list()
mod.20 <- list()
time <- list()

for(i in seq_along(outcome)){
  for(j in seq_along(inverse)){
    for(k in seq_along(lag)){
      if(j == 1){
        vac.var <- paste0("fullvax",lag[k])
      }
      if(j == 2){
        vac.var <- paste0("fullvax",lag[k],"i")
      }
      
      start.time <- Sys.time()
      res.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]] <- lmer(update(log(get(outcome[i])) ~ 1, 
                                                                                      paste0(" ~ . + ", vac.var, " + travel + pos_perc + SI + (1|PROVI)")),
                                                                               data = cov_prov)
      end.time <- Sys.time()
      time[["all"]][[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]] <- end.time - start.time
      
      mod.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]][["significant"]] <- ifelse(all(car::Anova(res.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]])[,3] < 0.05) == TRUE,
                                                                                                  1,0)
      mod.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]][["singular"]] <- isSingular(res.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]])
      mod.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]][["AIC"]] <- AIC(res.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]])
      
      start.time <- Sys.time()
      res.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]] <- lmer(update(log(get(outcome[i])) ~ 1, 
                                                                                     paste0(" ~ . + ", vac.var, " + travel + pos_perc + SI + (1|PROVI)")),
                                                                              data = cov_prov%>%filter(Week>=20))
      end.time <- Sys.time()
      time[["20"]][[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]] <- end.time - start.time
      
      mod.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]][["significant"]] <- ifelse(all(car::Anova(res.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]])[,3] < 0.05) == TRUE,
                                                                                                 1,0)
      mod.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]][["singular"]] <- isSingular(res.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]])
      mod.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]][["AIC"]] <- AIC(res.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]])
      
      for(l in c(1:6)){
        inter.term <- data.frame(X = c("vt", "vp", "vs", "tp", "ts", "ps"))%>%
          mutate(form = case_when(X == "vt" ~ paste0(vac.var,"*travel"),
                                  X == "vp" ~ paste0(vac.var,"*pos_perc"),
                                  X == "vs" ~ paste0(vac.var,"*SI"),
                                  X == "tp" ~ "travel*pos_perc",
                                  X == "ts" ~ "travel*SI",
                                  X == "ps" ~ "pos_perc*SI"))
        if(l == 1){
          fixef <- inter.term%>%
            rename(form.final = form)
        }
        if(l == 2){
          fixef <- data.frame(gtools::combinations(n = 6, r = 2, v = c("vt", "vp", "vs", "tp", "ts", "ps")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   X.final = paste0(X1,X2),
                   form.final = paste0(form.1,"+",form.2))
        }
        if(l == 3){
          fixef <- data.frame(gtools::combinations(n = 6, r = 3, v = c("vt", "vp", "vs", "tp", "ts", "ps")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   form.3 = inter.term$form[match(X3, inter.term$X)],
                   X.final = paste0(X1,X2,X3),
                   form.final = paste0(form.1,"+",form.2,"+",form.3))
        }
        if(l == 4){
          fixef <- data.frame(gtools::combinations(n = 6, r = 4, v = c("vt", "vp", "vs", "tp", "ts", "ps")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   form.3 = inter.term$form[match(X3, inter.term$X)],
                   form.4 = inter.term$form[match(X4, inter.term$X)],
                   X.final = paste0(X1,X2,X3,X4),
                   form.final = paste0(form.1,"+",form.2,"+",form.3,"+",form.4))
        }
        if(l == 5){
          fixef <- data.frame(gtools::combinations(n = 6, r = 5, v = c("vt", "vp", "vs", "tp", "ts", "ps")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   form.3 = inter.term$form[match(X3, inter.term$X)],
                   form.4 = inter.term$form[match(X4, inter.term$X)],
                   form.5 = inter.term$form[match(X5, inter.term$X)],
                   X.final = paste0(X1,X2,X3,X4,X5),
                   form.final = paste0(form.1,"+",form.2,"+",form.3,"+",form.4,"+",form.5))
        }
        if(l == 6){
          fixef <- data.frame(gtools::combinations(n = 6, r = 6, v = c("vt", "vp", "vs", "tp", "ts", "ps")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   form.3 = inter.term$form[match(X3, inter.term$X)],
                   form.4 = inter.term$form[match(X4, inter.term$X)],
                   form.5 = inter.term$form[match(X5, inter.term$X)],
                   form.6 = inter.term$form[match(X6, inter.term$X)],
                   X.final = paste0(X1,X2,X3,X4,X5,X6),
                   form.final = paste0(form.1,"+",form.2,"+",form.3,"+",form.4,"+",form.5,"+",form.6))
        }

        for(m in seq_along(fixef$X.final)){
          start.time <- Sys.time()
          res.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]] <- lmer(update(log(get(outcome[i])) ~ 1,
                                                                                                 paste0(" ~ . + ", vac.var, " + travel + pos_perc + SI + (1|PROVI) +",
                                                                                                        fixef$form.final[m])),
                                                                                          data = cov_prov)
          end.time <- Sys.time()
          time[["all"]][[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]] <- end.time - start.time

          mod.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]][["significant"]] <- ifelse(all(car::Anova(res.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]])[-c(1:4),3] < 0.05) == TRUE,
                                                                                                             1,0)
          mod.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]][["singular"]] <- isSingular(res.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]])
          mod.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]][["AIC"]] <- AIC(res.all[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]])

          start.time <- Sys.time()
          res.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]] <- lmer(update(log(get(outcome[i])) ~ 1,
                                                                                                paste0(" ~ . + ", vac.var, " + travel + pos_perc + SI + (1|PROVI) +",
                                                                                                       fixef$form.final[m])),
                                                                                         data = cov_prov%>%filter(Week>=20))
          end.time <- Sys.time()
          time[["20"]][[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]] <- end.time - start.time

          mod.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]][["significant"]] <- ifelse(all(car::Anova(res.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]])[-c(1:4),3] < 0.05) == TRUE,
                                                                                                            1,0)
          mod.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]][["singular"]] <- isSingular(res.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]])
          mod.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]][["AIC"]] <- AIC(res.20[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]])
        }
      }
    }
  }
}

res.all.alt <- list()
mod.all.alt <- list()
res.20.alt <- list()
mod.20.alt <- list()
time.alt <- list()

for(i in 2:2){
  for(j in seq_along(inverse)){
    for(k in seq_along(lag.alt)){
      if(j == 1){
        vac.var <- paste0("fullvax",lag.alt[k])
      }
      if(j == 2){
        vac.var <- paste0("fullvax",lag.alt[k],"i")
      }
      
      start.time <- Sys.time()
      res.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]] <- lmer(update(log(get(outcome[i])) ~ 1, 
                                                                                        paste0(" ~ . + ", vac.var, " + travel + pos_perc + SI + log_IC7 + (1|PROVI)")),
                                                                                 data = cov_prov)
      end.time <- Sys.time()
      time.alt[["all"]][[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]] <- end.time - start.time
      
      mod.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]][["significant"]] <- ifelse(all(car::Anova(res.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]])[,3] < 0.05) == TRUE,
                                                                                                    1,0)
      mod.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]][["singular"]] <- isSingular(res.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]])
      mod.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]][["AIC"]] <- AIC(res.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]])
      
      start.time <- Sys.time()
      res.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]] <- lmer(update(log(get(outcome[i])) ~ 1, 
                                                                                       paste0(" ~ . + ", vac.var, " + travel + pos_perc + SI + log_IC7 + (1|PROVI)")),
                                                                                data = cov_prov%>%filter(Week>=20))
      end.time <- Sys.time()
      time.alt[["20"]][[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]] <- end.time - start.time
      
      mod.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]][["significant"]] <- ifelse(all(car::Anova(res.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]])[,3] < 0.05) == TRUE,
                                                                                                   1,0)
      mod.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]][["singular"]] <- isSingular(res.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]])
      mod.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]][["AIC"]] <- AIC(res.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][["nointer"]])
      
      for(l in c(1:10)){
        inter.term <- data.frame(X = c("ip", "is", "it", "iv", "ps", "pt", "pv",
                                       "st", "sv", "tv"))%>%
          mutate(form = case_when(X == "ip" ~ "log_IC7*pos_perc",
                                  X == "is" ~ "log_IC7*SI",
                                  X == "it" ~ "log_IC7*travel",
                                  X == "iv" ~ paste0("log_IC7*",vac.var),
                                  X == "ps" ~ "pos_perc*SI",
                                  X == "pt" ~ "pos_perc*travel",
                                  X == "pv" ~ paste0("pos_perc*",vac.var),
                                  X == "st" ~ "SI*travel",
                                  X == "sv" ~ paste0("SI*",vac.var),
                                  X == "tv" ~ paste0("travel*",vac.var)))
        if(l == 1){
          fixef <- inter.term%>%
            rename(form.final = form)
        }
        if(l == 2){
          fixef <- data.frame(gtools::combinations(n = 10, r = 2, v = c("ip", "is", "it", "iv", "ps", "pt", "pv",
                                                                        "st", "sv", "tv")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   X.final = paste0(X1,X2),
                   form.final = paste0(form.1,"+",form.2))
        }
        if(l == 3){
          fixef <- data.frame(gtools::combinations(n = 10, r = 3, v = c("ip", "is", "it", "iv", "ps", "pt", "pv",
                                                                        "st", "sv", "tv")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   form.3 = inter.term$form[match(X3, inter.term$X)],
                   X.final = paste0(X1,X2,X3),
                   form.final = paste0(form.1,"+",form.2,"+",form.3))
        }
        if(l == 4){
          fixef <- data.frame(gtools::combinations(n = 10, r = 4, v = c("ip", "is", "it", "iv", "ps", "pt", "pv",
                                                                        "st", "sv", "tv")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   form.3 = inter.term$form[match(X3, inter.term$X)],
                   form.4 = inter.term$form[match(X4, inter.term$X)],
                   X.final = paste0(X1,X2,X3,X4),
                   form.final = paste0(form.1,"+",form.2,"+",form.3,"+",form.4))
        }
        if(l == 5){
          fixef <- data.frame(gtools::combinations(n = 10, r = 5, v = c("ip", "is", "it", "iv", "ps", "pt", "pv",
                                                                        "st", "sv", "tv")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   form.3 = inter.term$form[match(X3, inter.term$X)],
                   form.4 = inter.term$form[match(X4, inter.term$X)],
                   form.5 = inter.term$form[match(X5, inter.term$X)],
                   X.final = paste0(X1,X2,X3,X4,X5),
                   form.final = paste0(form.1,"+",form.2,"+",form.3,"+",form.4,"+",form.5))
        }
        if(l == 6){
          fixef <- data.frame(gtools::combinations(n = 10, r = 6, v = c("ip", "is", "it", "iv", "ps", "pt", "pv",
                                                                        "st", "sv", "tv")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   form.3 = inter.term$form[match(X3, inter.term$X)],
                   form.4 = inter.term$form[match(X4, inter.term$X)],
                   form.5 = inter.term$form[match(X5, inter.term$X)],
                   form.6 = inter.term$form[match(X6, inter.term$X)],
                   X.final = paste0(X1,X2,X3,X4,X5,X6),
                   form.final = paste0(form.1,"+",form.2,"+",form.3,"+",form.4,
                                       "+",form.5,"+",form.6))
        }
        if(l == 7){
          fixef <- data.frame(gtools::combinations(n = 10, r = 7, v = c("ip", "is", "it", "iv", "ps", "pt", "pv",
                                                                        "st", "sv", "tv")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   form.3 = inter.term$form[match(X3, inter.term$X)],
                   form.4 = inter.term$form[match(X4, inter.term$X)],
                   form.5 = inter.term$form[match(X5, inter.term$X)],
                   form.6 = inter.term$form[match(X6, inter.term$X)],
                   form.7 = inter.term$form[match(X7, inter.term$X)],
                   X.final = paste0(X1,X2,X3,X4,X5,X6,X7),
                   form.final = paste0(form.1,"+",form.2,"+",form.3,"+",form.4,
                                       "+",form.5,"+",form.6,"+",form.7))
        }
        if(l == 8){
          fixef <- data.frame(gtools::combinations(n = 10, r = 8, v = c("ip", "is", "it", "iv", "ps", "pt", "pv",
                                                                        "st", "sv", "tv")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   form.3 = inter.term$form[match(X3, inter.term$X)],
                   form.4 = inter.term$form[match(X4, inter.term$X)],
                   form.5 = inter.term$form[match(X5, inter.term$X)],
                   form.6 = inter.term$form[match(X6, inter.term$X)],
                   form.7 = inter.term$form[match(X7, inter.term$X)],
                   form.8 = inter.term$form[match(X8, inter.term$X)],
                   X.final = paste0(X1,X2,X3,X4,X5,X6,X7,X8),
                   form.final = paste0(form.1,"+",form.2,"+",form.3,"+",form.4,
                                       "+",form.5,"+",form.6,"+",form.7,"+",form.8))
        }
        if(l == 9){
          fixef <- data.frame(gtools::combinations(n = 10, r = 9, v = c("ip", "is", "it", "iv", "ps", "pt", "pv",
                                                                        "st", "sv", "tv")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   form.3 = inter.term$form[match(X3, inter.term$X)],
                   form.4 = inter.term$form[match(X4, inter.term$X)],
                   form.5 = inter.term$form[match(X5, inter.term$X)],
                   form.6 = inter.term$form[match(X6, inter.term$X)],
                   form.7 = inter.term$form[match(X7, inter.term$X)],
                   form.8 = inter.term$form[match(X8, inter.term$X)],
                   form.9 = inter.term$form[match(X9, inter.term$X)],
                   X.final = paste0(X1,X2,X3,X4,X5,X6,X7,X8,X9),
                   form.final = paste0(form.1,"+",form.2,"+",form.3,"+",form.4,
                                       "+",form.5,"+",form.6,"+",form.7,"+",form.8,
                                       "+",form.9))
        }
        if(l == 10){
          fixef <- data.frame(gtools::combinations(n = 10, r = 10, v = c("ip", "is", "it", "iv", "ps", "pt", "pv",
                                                                        "st", "sv", "tv")))%>%
            mutate(form.1 = inter.term$form[match(X1, inter.term$X)],
                   form.2 = inter.term$form[match(X2, inter.term$X)],
                   form.3 = inter.term$form[match(X3, inter.term$X)],
                   form.4 = inter.term$form[match(X4, inter.term$X)],
                   form.5 = inter.term$form[match(X5, inter.term$X)],
                   form.6 = inter.term$form[match(X6, inter.term$X)],
                   form.7 = inter.term$form[match(X7, inter.term$X)],
                   form.8 = inter.term$form[match(X8, inter.term$X)],
                   form.9 = inter.term$form[match(X9, inter.term$X)],
                   form.10 = inter.term$form[match(X10, inter.term$X)],
                   X.final = paste0(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10),
                   form.final = paste0(form.1,"+",form.2,"+",form.3,"+",form.4,
                                       "+",form.5,"+",form.6,"+",form.7,"+",form.8,
                                       "+",form.9,"+",form.10))
        }
        
        for(m in seq_along(fixef$X.final)){
          start.time <- Sys.time()
          res.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]] <- lmer(update(log(get(outcome[i])) ~ 1, 
                                                                                                   paste0(" ~ . + ", vac.var, " + travel + pos_perc + SI + (1|PROVI) +",
                                                                                                          fixef$form.final[m])),
                                                                                            data = cov_prov)
          end.time <- Sys.time()
          time.alt[["all"]][[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]] <- end.time - start.time
          
          mod.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]][["significant"]] <- ifelse(all(car::Anova(res.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]])[-c(1:4),3] < 0.05) == TRUE,
                                                                                                               1,0)
          mod.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]][["singular"]] <- isSingular(res.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]])
          mod.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]][["AIC"]] <- AIC(res.all.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]])
          
          start.time <- Sys.time()
          res.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]] <- lmer(update(log(get(outcome[i])) ~ 1, 
                                                                                                  paste0(" ~ . + ", vac.var, " + travel + pos_perc + SI + (1|PROVI) +",
                                                                                                         fixef$form.final[m])),
                                                                                           data = cov_prov%>%filter(Week>=20))
          end.time <- Sys.time()
          time.alt[["20"]][[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]] <- end.time - start.time
          
          mod.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]][["significant"]] <- ifelse(all(car::Anova(res.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]])[-c(1:4),3] < 0.05) == TRUE,
                                                                                                              1,0)
          mod.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]][["singular"]] <- isSingular(res.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]])
          mod.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]][["AIC"]] <- AIC(res.20.alt[[outcome[[i]]]][[inverse[j]]][[vac.var]][[fixef$X.final[m]]])
        }
      }
    }
  }
}
