### The variables in this analysis were structured as follows:
### ISOYEAR: yearly time interval based on ISO 8601 definition
### ISOWEEK: weekly time interval based on ISO 8601 definition
### CD_SEX: M = male, F = female
### CD_AGEGROUP: age categories (0-24, 25-44, 45-64, 65-74, 75-84, 85+)
### HEATWAVE: 0 = not part of heatwave, 1 = part of heatwave
### S1-S6: Fourier terms
### MR: mortality rate per 100,000


## Select reference data
mort.ref.young <- mort.all%>%
  filter(ISOYEAR >= 2009 & ISOYEAR < 2020 & CD_AGEGROUP %in% c("0-24", "25-44"))

mort.ref.old <- mort.all%>%
  filter(ISOYEAR >= 2009 & ISOYEAR < 2020 & !CD_AGEGROUP %in% c("0-24", "25-44"))

## Fit the reference models
lmm.2.2.young <- lmer(log(MR) ~ CD_AGEGROUP * CD_SEX +
                        (1|ISOYEAR), 
                      data = mort.ref.young, REML = TRUE)
summary(lmm.2.2.young)

lmm.2.2.old.H <- lmer(log(MR) ~ CD_AGEGROUP * CD_SEX + HEAT +
                        S1 + S2 + S3 + S4 + 
                        (1|ISOYEAR), 
                      data = mort.ref.old, REML = TRUE)
summary(lmm.2.2.old.H)

## Predict based on the fitted reference models
pred.data.young <- mort.all%>%
  filter(ISOYEAR %in% c(2020:2022) & CD_AGEGROUP %in% c("0-24", "25-44"))%>%
  as.data.frame()

pi.pred.young <- exp(predictInterval(merMod = lmm.2.2.young, newdata = pred.data.young,
                                     level = 0.95, n.sims = 1000,
                                     stat = "mean", type = "linear.prediction",
                                     include.resid.var = TRUE))

pred.data.old <- mort.all%>%
  filter(ISOYEAR %in% c(2020:2022) & !CD_AGEGROUP %in% c("0-24", "25-44"))%>%
  as.data.frame()

pi.pred.old <- exp(predictInterval(merMod = lmm.2.2.old.H, newdata = pred.data.old,
                                   level = 0.95, n.sims = 1000,
                                   stat = "mean", type = "linear.prediction",
                                   include.resid.var = TRUE))

## Plot the results together
fit.pred <- ggplot(data = bind_rows(bind_cols(pred.data.young, pi.pred.young),
                                    bind_cols(pred.data.old, pi.pred.old))%>%
                     mutate(y.upr = case_when(CD_AGEGROUP %in% c("0-24", "25-44") ~ 5,
                                              CD_AGEGROUP %in% c("45-64") ~ 20,
                                              CD_AGEGROUP %in% c("65-74") ~ 70,
                                              CD_AGEGROUP %in% c("75-84") ~ 250,
                                              CD_AGEGROUP %in% c("85+") ~ 800)), 
                   aes(x = ISOWEEK)) +
  geom_point(aes(y = MR, colour = as.factor(ISOYEAR)), 
             alpha = 0.8, size = 0.3) +
  geom_line(aes(y = fit), colour = "grey30", 
            linewidth = 0.3, linetype = "solid") +  # adding predicted values
  scale_colour_manual("Year", values = c("2020" = cbbPalette[10],
                                         "2021" = cbbPalette[11],
                                         "2022" = cbbPalette[4])) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey30", 
              alpha = 0.3, show.legend = FALSE) + # adding predicted intervals
  facet_wrap(factor(CD_SEX,
                    levels = c("M", "F"),
                    labels = c("Male", "Female")) ~ CD_AGEGROUP, 
             ncol = 6, scales = "free_y") +
  geom_blank(aes(y = 0)) +
  geom_blank(aes(y = y.upr)) +
  scale_x_continuous("Week", limits = c(1,52),
                     breaks = seq(1, 52, by = 10)) +
  scale_y_continuous("Mortality rate per 100,000 individuals") +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  theme_bw() +
  theme(plot.title = element_text(color = "black", size = 8, face = "bold"),
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.2, 'cm'),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        text = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_text(margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
                                  face = "bold"),
        axis.title.x = element_text(size = 8, color = "black", vjust = -2),
        axis.title.y = element_text(size = 8, color = "black", vjust = +2),
        axis.text = element_text(size = 8, color = "black"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))