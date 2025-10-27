### The variables in this analysis were structured as follows:
### YYYY: Year
### MM: Month
### CD_SEX: M = male, F = female
### CD_AGEGROUP: age categories (0-24, 25-44, 45-64, 65-74, 75-84, 85+)
### SI: median stringency index
### CHAP_COD_2: cause of death
### MR: mortality rate per 100,000

## Fit the model for each cause of death
COD_2 <- c("External causes", "Heart and vascular diseases", "Neoplasms", 
           "Other causes", "Infectious diseases", "Mental and behavioral disorders",
           "COVID-19")

res.cod <- list()
for(i in seq_along(COD_2)){
  res.cod[[COD_2[i]]][["1"]] <- lmer(log(MR) ~ CD_SEX * CD_AGEGROUP * SI + 
                                       (1|YYYY), 
                                     data = cod.all%>%
                                       filter(CHAP_COD_2 == COD_2[i]), REML = TRUE)
  
  res.cod[[COD_2[i]]][["2"]] <- lmer(log(MR) ~  CD_SEX*CD_AGEGROUP + CD_SEX*SI + 
                                       CD_AGEGROUP*SI +
                                       (1|YYYY), 
                                     data = cod.all%>%
                                       filter(CHAP_COD_2 == COD_2[i]), REML = TRUE)
  
  res.cod[[COD_2[i]]][["3"]] <- lmer(log(MR) ~  CD_SEX*CD_AGEGROUP + CD_SEX*SI + 
                                       (1|YYYY), 
                                     data = cod.all%>%
                                       filter(CHAP_COD_2 == COD_2[i]), REML = TRUE)
  
  res.cod[[COD_2[i]]][["4"]] <- lmer(log(MR) ~  CD_SEX*CD_AGEGROUP + CD_AGEGROUP*SI +
                                       (1|YYYY), 
                                     data = cod.all%>%
                                       filter(CHAP_COD_2 == COD_2[i]), REML = TRUE)
  
  res.cod[[COD_2[i]]][["5"]] <- lmer(log(MR) ~  CD_SEX*SI + CD_AGEGROUP*SI +
                                       (1|YYYY), 
                                     data = cod.all%>%
                                       filter(CHAP_COD_2 == COD_2[i]), REML = TRUE)
  
  res.cod[[COD_2[i]]][["6"]] <- lmer(log(MR) ~ CD_SEX * CD_AGEGROUP + SI + 
                                       (1|YYYY), 
                                     data = cod.all%>%
                                       filter(CHAP_COD_2 == COD_2[i]), REML = TRUE)
  
  res.cod[[COD_2[i]]][["7"]] <- lmer(log(MR) ~ CD_SEX * SI + CD_AGEGROUP + 
                                       (1|YYYY), 
                                     data = cod.all%>%
                                       filter(CHAP_COD_2 == COD_2[i]), REML = TRUE)
  
  res.cod[[COD_2[i]]][["8"]] <- lmer(log(MR) ~ CD_SEX + CD_AGEGROUP * SI + 
                                       (1|YYYY), 
                                     data = cod.all%>%
                                       filter(CHAP_COD_2 == COD_2[i]), REML = TRUE)
  
  res.cod[[COD_2[i]]][["9"]] <- lmer(log(MR) ~ CD_SEX + CD_AGEGROUP + SI + 
                                       (1|YYYY), 
                                     data = cod.all%>%
                                       filter(CHAP_COD_2 == COD_2[i]), REML = TRUE)
}