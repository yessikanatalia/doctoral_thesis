rm(list = ls())
pacman::p_load(dplyr, tidyr, fractaldim, zoo, forecast,
               cowplot, ggplot2, scales, ggforce,
               factoextra)

### Function for fractal dimension
fract.d = function (data, window, method) {
  rollapply(data, width = window, 
            FUN = function (x) {
              return(tryCatch(fd.estimate(x, method = method, keep.loglog = T,
                                          plot.loglog = F, plot.allpoints = T,
                                          nlags = "auto")$fd, 
                              error=function(e) NA))},
            align = "right", partial = T)
  
}

fillNA = function(x) {
  ifelse(is.na(x),1,x)
}

### Wrapper function for replications
cases_sim <- function (N, Lambda, Methods, Window){
  nbase <- list()
  for (i in seq_along(Lambda)){
    data0 <- data.frame(cbind(Day = c(1:N), 
                              n0 = arima.sim(list(order = c(0,0,0)), 
                                             rand.gen = function(x) rpois(x, lambda = Lambda[i]),
                                             n = N)))
    nbase[[i]] <- data.frame(crossing(Day = c(1:N), Methods, Window, Lambda[i]))%>%
      left_join(., data0, by = "Day")%>%
      group_by(Methods, Window)%>%
      arrange(Day)%>%
      mutate(fd = fillNA(fract.d(n0, window = as.numeric(unique(Window)), 
                                 method = unique(Methods))))
  }
  nbase_long <- bind_rows(nbase)%>%
    rename(Lambda = Lambda.i.)
  
  sumbase <- nbase_long%>%
    group_by(Methods, Window, Lambda)%>%
    summarise(mean.fd = round(mean(fd, na.rm = T),2),
              variance.fd = round(var(fd),2),
              acf.fd = round(acf(fd,lag.max = 1, plot = F, na.action = na.pass)$acf[2],2))
  sumbase
}

### Generate one example 
lamb_pois <- seq(0.01,0.5,0.05) #lambda for the Poisson distribution
n <- 365 #length of series
Day <- c(1:n)
Methods <- c("boxcount", "hallwood", "variogram", "madogram")
Window <- c(7,14,21)

n_list <- list()
p_list <- list()
for (i in seq_along(lamb_pois)){
  set.seed(0)
  nbase <- data.frame(cbind(Day, n0 = arima.sim(list(order = c(0,0,0)), 
                                                rand.gen = function(x) rpois(x, lambda = lamb_pois[i]),
                                                n = n)))
  n_list[[i]] <- data.frame(crossing(Day, Methods, Window))%>%
    left_join(., nbase, by = "Day")%>%
    mutate(sim = i)%>%
    group_by(Methods, Window)%>%
    arrange(Day)%>%
    mutate(fd = fillNA(fract.d(n0, window = as.numeric(unique(Window)), 
                               method = unique(Methods))))
  
  p.fd <- ggplot(n_list[[i]], aes(x = Day)) + 
    geom_smooth(aes(y = fd, color = Methods), method = "loess", 
                span = 0.05, linewidth = 0.3, se = FALSE, method.args = list(degree = 1)) +
    scale_color_manual(values = c("boxcount" = "#CC6677",
                                  "hallwood" = "#117733",
                                  "variogram" = "#DDCC77",
                                  "madogram" = "#88CCEE"),
                       labels = c("boxcount" = "Boxcount",
                                  "hallwood" = "Hall-Wood",
                                  "variogram" = "Variogram",
                                  "madogram" = "Madogram")) +
    scale_x_continuous(name = "") +
    scale_y_continuous("Fractal dimension", limits = c(1,2)) +
    facet_wrap(~ factor(Window, levels = c("7", "14", "21")),
               ncol = 1) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 7),
          legend.background = element_rect(fill = NA),
          legend.box.background = element_rect(fill = NA, colour = NA),
          axis.text.y = element_text(size = 7, color = c("black")),
          axis.ticks.y = element_line(color = c("black")),
          axis.title = element_text(size = 7, color = c("black")),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          text = element_text(size = 7),
          strip.text = element_text(margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")),
          plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"))
  
  p.legend <- cowplot::get_plot_component(p.fd, 'guide-box', return_all = TRUE)[[3]]
  
  p.fd.final <- p.fd + theme(legend.position = "none")
  
  p.cases <- ggplot(nbase, aes(x = Day)) + 
    geom_line(aes(y = n0), linewidth = 0.3, color = "black")+
    scale_x_continuous(name = "Day") +
    scale_y_continuous("Incidence", limits = c(0,4)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 7, color = c("black")),
          axis.ticks.y = element_line(color = c("black")),
          axis.title = element_text(size = 7, color = c("black")),
          axis.text.x = element_text(size = 7, color = "black"),
          text = element_text(size = 7),
          plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")) 
  
  title <- ggdraw() + 
    draw_label(paste("\u03BB","=",lamb_pois[i]), fontface = 'bold', size = 7) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 
  
  p_list[[i]] <- cowplot::plot_grid(title, NULL, p.fd.final, p.cases,
                                    ncol = 1, align = "v", axis = "tblr",
                                    rel_heights = c(0.2,-0.05,1,0.4))
}

ggsave(cowplot::plot_grid(cowplot::plot_grid(p_list[[1]], p_list[[2]], 
                                             p_list[[3]], p_list[[4]],
                                             ncol = 2), 
                          p.legend, ncol = 1, rel_heights = c(2,0.1)),
       file = paste0("G:/My Drive/COVID-19-Fractal dimension/figures/C1.Sim1.png"), 
       width = 2000, height = 2000, units = "px", device = "png", dpi = 300)

ggsave(cowplot::plot_grid(cowplot::plot_grid(p_list[[5]], p_list[[6]], 
                                             p_list[[7]], p_list[[8]],
                                             p_list[[9]], p_list[[10]],
                                             ncol = 2),
                          p.legend, ncol = 1, rel_heights = c(2,0.1)),
       file = paste0("G:/My Drive/COVID-19-Fractal dimension/figures/C1.Sim2.png"), 
       width = 2000, height = 3000, units = "px", device = "png", dpi = 300)

### Generate 1000 replications
N <- 365
lamb_pois <- seq(0.01,0.5,0.05)
methods <- c("boxcount", "hallwood", "variogram", "madogram")
window <- c(7,14,21)

sim_list <- list()
set.seed(0)
for (i in 1:1000){
  sim_list[[i]] <- cases_sim(N = N, Lambda = lamb_pois, Methods = methods, Window = window)
  sim_list[[i]]$sim <- i
}

sim_long <- bind_rows(sim_list)%>%
  mutate(meth.class = case_when(Methods== "boxcount" ~ "Boxcount",
                                Methods== "hallwood" ~ "Hall-Wood",
                                Methods== "variogram" ~ "Variogram",
                                Methods== "madogram" ~ "Madogram"))

box_list <- list()
# lamb_select <- c(lamb_pois[1], lamb_pois[5], lamb_pois[10]) # selected lambda's
for (i in seq_along(lamb_pois))
{
  data0 <- sim_long%>%
    filter(Lambda == lamb_pois[i])
  
  p.mean <- ggplot(data0) +
    geom_boxplot(aes(x = meth.class, y = mean.fd, fill = meth.class),
                 width = 0.3, linewidth = 0.3, size = 0.3, outlier.size = 0.1) +
    scale_fill_manual(values = c("Boxcount" = "#CC6677",
                                 "Hall-Wood" = "#117733",
                                 "Variogram" = "#DDCC77",
                                 "Madogram" = "#88CCEE")) +
    scale_y_continuous("Mean FD", limits = c(1,2)) +
    labs(x = "") +
    facet_wrap(~ factor(Window, levels = c("7", "14", "21")),
               ncol = 1) +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 7),
          axis.text.x = element_text(size = 7,colour = "black"),
          axis.text.y.left = element_text(size = 7,colour = "black"),
          strip.text = element_text(margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")),
          plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm"))
  
  p.var <- ggplot(data0) +
    geom_boxplot(aes(x = meth.class, y = variance.fd, fill = meth.class),
                 width = 0.3, linewidth = 0.3, size = 0.3, outlier.size = 0.1) +
    scale_fill_manual(values = c("Boxcount" = "#CC6677",
                                 "Hall-Wood" = "#117733",
                                 "Variogram" = "#DDCC77",
                                 "Madogram" = "#88CCEE")) +
    scale_y_continuous("Variance FD", limits = c(0,0.3)) +
    labs(x = "") +
    facet_wrap(~ factor(Window, levels = c("7", "14", "21")),
               ncol = 1) +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 7),
          axis.text.x = element_text(size = 7,colour = "black"),
          axis.text.y.left = element_text(size = 7,colour = "black"),
          strip.text = element_text(margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")),
          plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm"))
  
  p.acf <- ggplot(data0) +
    geom_boxplot(aes(x = meth.class, y = acf.fd, fill = meth.class),
                 width = 0.3, linewidth = 0.3, size = 0.3, outlier.size = 0.1) +
    scale_fill_manual(values = c("Boxcount" = "#CC6677",
                                 "Hall-Wood" = "#117733",
                                 "Variogram" = "#DDCC77",
                                 "Madogram" = "#88CCEE")) +
    scale_y_continuous("ACF FD", limits = c(-1,1)) +
    labs(x = "") +
    facet_wrap(~ factor(Window, levels = c("7", "14", "21")),
               ncol = 1) +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 7),
          axis.text.x = element_text(size = 7,colour = "black"),
          axis.text.y.left = element_text(size = 7,colour = "black"),
          strip.text = element_text(margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")),
          plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm"))
  
  title <- ggdraw() + 
    draw_label(paste("\u03BB","=",lamb_pois[i]), fontface = "bold", size = 7)
  
  box_list[[i]] <- cowplot::plot_grid(title, NULL, p.mean,p.var,p.acf,
                                      ncol = 1, align = "v", axis = "tblr",
                                      rel_heights = c(0.2,-0.05,1.1,1.1,1.1))
}

# Create individual figures
ggsave(cowplot::plot_grid(box_list[[1]],box_list[[2]],box_list[[3]],box_list[[4]],box_list[[5]],
                          ncol = 5),
       file = paste0("G:/My Drive/COVID-19-Fractal dimension/figures/C2.Sim1.png"), 
       width = 4000, height = 2500, units = "px", device = "png", dpi = 300)

ggsave(cowplot::plot_grid(box_list[[6]],box_list[[7]],box_list[[8]],box_list[[9]],box_list[[10]],
                          ncol = 5),
       file = paste0("G:/My Drive/COVID-19-Fractal dimension/figures/C2.Sim2.png"), 
       width = 4000, height = 2500, units = "px", device = "png", dpi = 300)

### Cluster analysis (k-means)
## Check for optimal number of cluster
for (j in seq_along(methods)){
  for(k in seq_along(window)){
    subset <- sim_long%>%
      filter(sim %in% c(1:5) & Methods == methods[j] & Window == window[k])%>%
      ungroup()%>%
      dplyr::select(mean.fd, variance.fd, acf.fd)%>%
      na.omit()%>%
      data.frame()
    
    print(fviz_nbclust(subset, kmeans, method = "wss") +
            labs(subtitle = paste0(methods[j], ", sliding window = ", window[k], " days")))
  }
}

## Create cluster figures
data.clust <- list()
fig.clust <- list()

# All replications together
for (j in seq_along(methods)){
    for(k in seq_along(window)){
      subset.1 <- sim_long%>%
        filter(sim %in% c (1:10) & Methods == methods[j] & Window == window[k])%>%
        ungroup()
      
      subset.2 <- subset.1%>%
        dplyr::select(mean.fd, variance.fd, acf.fd)%>%
        na.omit()%>%
        data.frame()
      
      # Calculate the clusters
      clust.0 <- kmeans(subset.2, centers = 4, nstart = 25)
      data.clust[[methods[j]]][[paste(window[k])]] <- bind_cols(subset.1, 
                                                                            cluster = as.character(clust.0[["cluster"]]))%>%
        left_join(., tibble::rownames_to_column(data.frame(clust.0[["centers"]]), "cluster")%>%
                    rename(mean.centroid = mean.fd,
                           variance.centroid = variance.fd,
                           acf.centroid = acf.fd),
                  by = "cluster")%>%
        group_by(Methods, Window)%>%
        mutate(mean.mean = mean(mean.fd, na.rm = T),
               mean.variance = mean(variance.fd, na.rm = T),
               mean.acf = mean(acf.fd, na.rm = T))%>%
        ungroup()%>%
        mutate(clust.labs = case_when(mean.centroid < mean.mean & variance.centroid < mean.variance &
                                        acf.centroid < mean.acf ~ "Low mean, low variance, low ACF",
                                      mean.centroid < mean.mean &  variance.centroid < mean.variance &
                                        acf.centroid >= mean.acf ~ "Low mean, low variance, high ACF",
                                      mean.centroid < mean.mean & variance.centroid >= mean.variance & 
                                        acf.centroid < mean.acf ~ "Low mean, high variance, low ACF",
                                      mean.centroid < mean.mean & variance.centroid >= mean.variance &
                                        acf.centroid >= mean.acf ~ "Low mean, high variance, high ACF",
                                      mean.centroid >= mean.mean & variance.centroid < mean.variance &
                                        acf.centroid < mean.acf ~ "High mean, low variance, low ACF",
                                      mean.centroid >= mean.mean & variance.centroid < mean.variance &
                                        acf.centroid >= mean.acf ~ "High mean, low variance, high ACF",
                                      mean.centroid >= mean.mean & variance.centroid >= mean.variance &
                                        acf.centroid < mean.acf ~ "High mean, high variance, low ACF",
                                      mean.centroid >= mean.mean & variance.centroid >= mean.variance & 
                                        acf.centroid >= mean.acf ~ "High mean, high variance, high ACF"))
      
      # Create initial figures
      fig.clust[[methods[j]]][[paste(window[k])]][["init"]] <- fviz_cluster(clust.0, data = subset.2, 
                                                                                        geom = "point", pointsize = 0.2) +
        scale_color_brewer("", palette = "Set1",
                           labels = c("1" = unique(data.clust[[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[methods[j]]][[paste(window[k])]]$cluster == "1"]),
                                      "2" = unique(data.clust[[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[methods[j]]][[paste(window[k])]]$cluster == "2"]),
                                      "3" = unique(data.clust[[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[methods[j]]][[paste(window[k])]]$cluster == "3"]),
                                      "4" = unique(data.clust[[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[methods[j]]][[paste(window[k])]]$cluster == "4"]))) + 
        scale_fill_brewer("", palette = "Set1",
                          labels = c("1" = unique(data.clust[[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[methods[j]]][[paste(window[k])]]$cluster == "1"]),
                                     "2" = unique(data.clust[[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[methods[j]]][[paste(window[k])]]$cluster == "2"]),
                                     "3" = unique(data.clust[[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[methods[j]]][[paste(window[k])]]$cluster == "3"]),
                                     "4" = unique(data.clust[[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[methods[j]]][[paste(window[k])]]$cluster == "4"]))) +
        scale_shape_manual("", values = c(22,23,24,25), 
                           labels = c("1" = unique(data.clust[[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[methods[j]]][[paste(window[k])]]$cluster == "1"]),
                                      "2" = unique(data.clust[[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[methods[j]]][[paste(window[k])]]$cluster == "2"]),
                                      "3" = unique(data.clust[[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[methods[j]]][[paste(window[k])]]$cluster == "3"]),
                                      "4" = unique(data.clust[[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[methods[j]]][[paste(window[k])]]$cluster == "4"]))) + 
        labs(title = paste0(window[k])) +
        guides(fill = guide_legend(nrow = 4), color = guide_legend(nrow = 4), 
               shape = guide_legend(nrow = 4)) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "bottom",
              legend.text = element_text(size = 7),
              legend.margin = margin(r = 0,l = 0, t = 0, b = 0),
              legend.key.size = unit(0.25, 'cm'),
              legend.spacing.y = unit(0.1, 'cm'),
              legend.key.spacing.y = unit(0.1, 'cm'),
              text = element_text(size = 7),
              axis.text.x = element_text(size = 7,colour = "black"),
              axis.text.y.left = element_text(size = 7,colour = "black"),
              strip.text = element_text(margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")),
              plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
      
      fig.clust[[methods[j]]][[paste(window[k])]][["final"]] <- fig.clust[[methods[j]]][[paste(window[k])]][["init"]] + 
        geom_text(data = bind_cols(fig.clust[[methods[j]]][[paste(window[k])]][["init"]]$data,
                                   subset.1), 
                  aes(x = x, y = y, label = Lambda, color = cluster),
                  size = 1.5, vjust = -1, show.legend = F)
    }
    
    title <- ggdraw() + 
      draw_label(unique(subset.1$meth.class), fontface = 'bold', size = 7) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    
    fig.clust[[methods[j]]][["final"]] <- cowplot::plot_grid(title, NULL, 
                                                                         fig.clust[[methods[j]]][[paste(window[1])]][["final"]],
                                                                         fig.clust[[methods[j]]][[paste(window[2])]][["final"]],
                                                                         fig.clust[[methods[j]]][[paste(window[3])]][["final"]],
                                                                         ncol = 1, align = "v", axis = "tblr",
                                                                         rel_heights = c(0.2,-0.05,1.1,1.1,1.1))
}

ggsave(cowplot::plot_grid(fig.clust[[methods[1]]][["final"]],
                          fig.clust[[methods[2]]][["final"]],
                          fig.clust[[methods[3]]][["final"]],
                          fig.clust[[methods[4]]][["final"]],
                          ncol = 4, align = "v", axis = "tblr"),
       file = paste0("G:/My Drive/COVID-19-Fractal dimension/figures/C3.all.png"),
       width = 2700, height = 2700, units = "px", device = 'png', dpi = 300)

# Per replication
for(i in 1:10){
  for (j in seq_along(methods)){
    for(k in seq_along(window)){
      subset.1 <- sim_long%>%
        filter(sim == i & Methods == methods[j] & Window == window[k])%>%
        ungroup()
      
      subset.2 <- subset.1%>%
        dplyr::select(mean.fd, variance.fd, acf.fd)%>%
        na.omit()%>%
        data.frame()
      
      # Calculate the clusters
      clust.0 <- kmeans(subset.2, centers = 4, nstart = 25)
      data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]] <- bind_cols(subset.1, 
                                                                            cluster = as.character(clust.0[["cluster"]]))%>%
        left_join(., tibble::rownames_to_column(data.frame(clust.0[["centers"]]), "cluster")%>%
                    rename(mean.centroid = mean.fd,
                           variance.centroid = variance.fd,
                           acf.centroid = acf.fd),
                  by = "cluster")%>%
        group_by(Methods, Window)%>%
        mutate(mean.mean = mean(mean.fd, na.rm = T),
               mean.variance = mean(variance.fd, na.rm = T),
               mean.acf = mean(acf.fd, na.rm = T))%>%
        ungroup()%>%
        mutate(clust.labs = case_when(mean.centroid < mean.mean & variance.centroid < mean.variance &
                                        acf.centroid < mean.acf ~ "Low mean, low variance, low ACF",
                                      mean.centroid < mean.mean &  variance.centroid < mean.variance &
                                        acf.centroid >= mean.acf ~ "Low mean, low variance, high ACF",
                                      mean.centroid < mean.mean & variance.centroid >= mean.variance & 
                                        acf.centroid < mean.acf ~ "Low mean, high variance, low ACF",
                                      mean.centroid < mean.mean & variance.centroid >= mean.variance &
                                        acf.centroid >= mean.acf ~ "Low mean, high variance, high ACF",
                                      mean.centroid >= mean.mean & variance.centroid < mean.variance &
                                        acf.centroid < mean.acf ~ "High mean, low variance, low ACF",
                                      mean.centroid >= mean.mean & variance.centroid < mean.variance &
                                        acf.centroid >= mean.acf ~ "High mean, low variance, high ACF",
                                      mean.centroid >= mean.mean & variance.centroid >= mean.variance &
                                        acf.centroid < mean.acf ~ "High mean, high variance, low ACF",
                                      mean.centroid >= mean.mean & variance.centroid >= mean.variance & 
                                        acf.centroid >= mean.acf ~ "High mean, high variance, high ACF"))
      
      # Create initial figures
      fig.clust[[paste(i)]][[methods[j]]][[paste(window[k])]][["init"]] <- fviz_cluster(clust.0, data = subset.2, 
                                                                                        geom = "point", pointsize = 0.2) +
        scale_color_brewer("", palette = "Set1",
                           labels = c("1" = unique(data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$cluster == "1"]),
                                      "2" = unique(data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$cluster == "2"]),
                                      "3" = unique(data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$cluster == "3"]),
                                      "4" = unique(data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$cluster == "4"]))) + 
        scale_fill_brewer("", palette = "Set1",
                          labels = c("1" = unique(data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$cluster == "1"]),
                                     "2" = unique(data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$cluster == "2"]),
                                     "3" = unique(data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$cluster == "3"]),
                                     "4" = unique(data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$cluster == "4"]))) +
        scale_shape_manual("", values = c(22,23,24,25), 
                           labels = c("1" = unique(data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$cluster == "1"]),
                                      "2" = unique(data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$cluster == "2"]),
                                      "3" = unique(data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$cluster == "3"]),
                                      "4" = unique(data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$clust.labs[data.clust[[paste(i)]][[methods[j]]][[paste(window[k])]]$cluster == "4"]))) + 
        labs(title = paste0(window[k])) +
        guides(fill = guide_legend(nrow = 4),
               color = guide_legend(nrow = 4), 
               shape = guide_legend(nrow = 4)) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "bottom",
              legend.text = element_text(size = 7),
              legend.margin = margin(r = 0,l = 0, t = 0, b = 0),
              legend.key.size = unit(0.25, 'cm'),
              text = element_text(size = 7),
              axis.text.x = element_text(size = 7,colour = "black"),
              axis.text.y.left = element_text(size = 7,colour = "black"),
              strip.text = element_text(margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")),
              plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
      
      fig.clust[[paste(i)]][[methods[j]]][[paste(window[k])]][["final"]] <- fig.clust[[paste(i)]][[methods[j]]][[paste(window[k])]][["init"]] + 
        geom_text(data = bind_cols(fig.clust[[paste(i)]][[methods[j]]][[paste(window[k])]][["init"]]$data,
                                   subset.1), 
                  aes(x = x, y = y, label = Lambda, color = cluster),
                  size = 2, vjust = -1, show.legend = F)
    }
    
    title <- ggdraw() + 
      draw_label(unique(subset.1$meth.class), fontface = 'bold', size = 7) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    
    fig.clust[[paste(i)]][[methods[j]]][["final"]] <- cowplot::plot_grid(title, NULL, 
                                                                         fig.clust[[paste(i)]][[methods[j]]][[paste(window[1])]][["final"]],
                                                                         fig.clust[[paste(i)]][[methods[j]]][[paste(window[2])]][["final"]],
                                                                         fig.clust[[paste(i)]][[methods[j]]][[paste(window[3])]][["final"]],
                                                                         ncol = 1, align = "v", axis = "tblr",
                                                                         rel_heights = c(0.2,-0.05,1.1,1.1,1.1))
  }
}

ggsave(cowplot::plot_grid(fig.clust[[paste(1)]][[methods[1]]][["final"]],
                          fig.clust[[paste(1)]][[methods[2]]][["final"]],
                          fig.clust[[paste(1)]][[methods[3]]][["final"]],
                          fig.clust[[paste(1)]][[methods[4]]][["final"]],
                          ncol = 4, align = "v", axis = "tblr"),
       file = paste0("G:/My Drive/COVID-19-Fractal dimension/figures/C3.",1,".png"),
       width = 3000, height = 3000, units="px", device='png', dpi=300)
