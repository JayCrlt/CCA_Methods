#### Set up environment ----

## Global options
rm(list = ls()) ; options(warn = -1, cores = parallel::detectCores()) ; setwd("..")

## Packages
library(geomtextpath)
library(viridisLite)
library(hrbrthemes)
library(tidybayes)
library(patchwork)
library(posterior)
library(tidyverse) 
library(pkgconfig)
library(ggridges) 
library(viridis)
library(ggimage)
library(Matrix)
library(readxl)
library(brms)
library(rsvg)

## Loading systematic review analyse
source("R/Scripts/Systematic_Review_CCA.R")

## Loading images
massive_svg   = data.frame(x = 80, y = 50, image = "Figures/Icons/platygyra-spp-brain-coral.svg")
branching_svg = data.frame(x = 80, y = 50, image = "Figures/Icons/acropora-cervicornis-staghorn-coral.svg")
coral_cca_svg = data.frame(x = 80, y = 50, image = "Figures/Icons/plate-coral-cca.svg")

## Loading Models
load("Models/fit_extension.RData")
load("Models/fit_surface.RData")
load("Models/fit_biom.RData")

## Empty vectors (to store object)
GR_viz        = vector("list", 8)
vector_method = vector("list", 8)
summary_list  = vector("list", 3)
color_list    = vector("list", 3)
rugosity.list = vector("list", 400)

## Color gradient
colors                           <- c("#ff595e", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93")
rates_plot_manual_legend_colours <- c("#a1dab4", "#fee090", "#f46d43", "#41b6c4")

# Useful functions
'%notin%' <- function(x,y)!('%in%'(x,y))

#### Figure 2 – Calcification rates CCA vs. corals ----

# Changing unit
data_cca = Data_viz[[3]] %>% mutate(., Rate_std = Rate_std/1000*365.25, std_error = std_error/1000*365.25)

# Building aan average dataset + a species dataset w/ quantiles informations
## CCA
Quantile_CCA_Data_all = data_cca %>% dplyr::filter(Climate %in% c("Warm temperate", "Tropical")) %>% 
  group_by(Genus) %>% summarise(.,         Q_0.01 = quantile(Rate_std, probs = 0.01),
                                Q_0.25 = quantile(Rate_std, probs = 0.25),
                                Q_0.50 = quantile(Rate_std, probs = 0.50),
                                Q_0.75 = quantile(Rate_std, probs = 0.75),
                                Q_0.99 = quantile(Rate_std, probs = 0.99))

Quantile_CCA_Data_avg = data.frame(Genus = "All CCA", 
                                   data_cca %>% dplyr::filter(Climate %in% c("Warm temperate", "Tropical")) %>% 
                                     dplyr::filter(Genus %notin% data_articulate$Genus) %>% 
                                     dplyr::filter(Genus != "Halimeda") %>% 
                                     group_by(Genus) %>% 
                                     summarise(., Q_0.01 = quantile(Rate_std, probs = 0.01),
                                               Q_0.25 = quantile(Rate_std, probs = 0.25),
                                               Q_0.50 = quantile(Rate_std, probs = 0.50),
                                               Q_0.75 = quantile(Rate_std, probs = 0.75),
                                               Q_0.99 = quantile(Rate_std, probs = 0.99)) %>% 
                                     summarise(., Q_0.01 = quantile(Q_0.01, probs = 0.01),
                                               Q_0.25 = quantile(Q_0.25, probs = 0.25),
                                               Q_0.50 = quantile(Q_0.50, probs = 0.50),
                                               Q_0.75 = quantile(Q_0.75, probs = 0.75),
                                               Q_0.99 = quantile(Q_0.99, probs = 0.99)))

## Corals
Quantile_Cor_Data_all = data_Niklas %>% group_by(Genus) %>% summarise(., 
                               Q_0.01 = quantile(conX, probs = 0.01),
                               Q_0.25 = quantile(conX, probs = 0.25),
                               Q_0.50 = quantile(conX, probs = 0.50),
                               Q_0.75 = quantile(conX, probs = 0.75),
                               Q_0.99 = quantile(conX, probs = 0.99))

Quantile_Cor_Data_avg = data.frame(Genus = "All Corals", data_Niklas %>% 
                                     group_by(Genus) %>% 
                                     summarise(., Q_0.01 = quantile(conX, probs = 0.01),
                                               Q_0.25 = quantile(conX, probs = 0.25),
                                               Q_0.50 = quantile(conX, probs = 0.50),
                                               Q_0.75 = quantile(conX, probs = 0.75),
                                               Q_0.99 = quantile(conX, probs = 0.99)) %>% 
                                     summarise(., Q_0.01 = quantile(Q_0.01, probs = 0.01),
                                               Q_0.25 = quantile(Q_0.25, probs = 0.25),
                                               Q_0.50 = quantile(Q_0.50, probs = 0.50),
                                               Q_0.75 = quantile(Q_0.75, probs = 0.75),
                                               Q_0.99 = quantile(Q_0.99, probs = 0.99)))

# Combine into 2 unique datasets
quantile_CCA_all = rbind(Quantile_CCA_Data_all, Quantile_CCA_Data_avg)
Quantile_Cor_all = rbind(Quantile_Cor_Data_all, Quantile_Cor_Data_avg)

# Combine into a single unique dataset
summary_all <- rbind(Quantile_Cor_all %>% dplyr::filter(., Genus == "All Corals"),
                     quantile_CCA_all %>% dplyr::filter(., Genus == "All CCA"))

# Starting to plot Figure 2 – Calcification rate CCA vs corals

Figure_2A <- summary_all %>%
  ggplot(aes(x = Genus, y = Q_0.50, col = Genus), show.legend = F) + 
  geom_linerange(aes(ymin = Q_0.01, ymax = Q_0.99), size = 1, show.legend = F) +
  geom_linerange(aes(ymin = Q_0.25, ymax = Q_0.75), size = 3, show.legend = F) + theme_bw() +
  geom_point(size = 5, show.legend = F, aes(shape = Genus, fill = Genus), col = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_shape_manual(name = "", limits = c("All Corals", "All CCA"), values = c(23, 21)) +
  scale_color_manual(values = c("violetred", "orange")) + scale_fill_manual(values = c("violetred", "orange")) +
  scale_y_continuous(name = expression("Calcification rate (g."*cm^-2*".yr"^-1*")"), 
                     limits = c(0,1), breaks = seq(0,1,0.2))

Figure_2Ba <- quantile_CCA_all %>% dplyr::filter(., Genus != "All CCA", Genus %notin% c("Amphiroa", "Halimeda")) %>% 
  ggplot(aes(x = Genus, y = Q_0.50)) + 
  geom_linerange(aes(ymin = Q_0.01, ymax = Q_0.99), size = 1, col = "violetred") +
  geom_linerange(aes(ymin = Q_0.25, ymax = Q_0.75), size = 3, col = "violetred") + theme_bw() +
  geom_point(size = 5, show.legend = F, shape = 21, fill = "violetred") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_x_discrete(name = "") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "", limits = c(0,1), breaks = seq(0,1,0.1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

Figure_2Bb <- quantile_CCA_all %>% dplyr::filter(., Genus == "Amphiroa") %>% 
  ggplot(aes(x = Genus, y = Q_0.50)) + 
  geom_linerange(aes(ymin = Q_0.01, ymax = Q_0.99), size = 1, col = "darkred") +
  geom_linerange(aes(ymin = Q_0.25, ymax = Q_0.75), size = 3, col = "darkred") + theme_bw() +
  geom_point(size = 5, show.legend = F, shape = 21, fill = "darkred") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_x_discrete(name = "") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "", limits = c(0,1), breaks = seq(0,1,0.1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

Figure_2Bc <- quantile_CCA_all %>% dplyr::filter(., Genus == "Halimeda") %>% 
  ggplot(aes(x = Genus, y = Q_0.50)) + 
  geom_linerange(aes(ymin = Q_0.01, ymax = Q_0.99), size = 1, col = "darkolivegreen4") +
  geom_linerange(aes(ymin = Q_0.25, ymax = Q_0.75), size = 3, col = "darkolivegreen4") + theme_bw() +
  geom_point(size = 5, show.legend = F, shape = 21, fill = "darkolivegreen4") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_x_discrete(name = "") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "", limits = c(0,1), breaks = seq(0,1,0.1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

Figure_2C <- Quantile_Cor_all %>% dplyr::filter(., Genus != "All Corals") %>% 
  ggplot(aes(x = Genus, y = Q_0.50)) + 
  geom_linerange(aes(ymin = Q_0.01, ymax = Q_0.99), size = 1, col = "orange") +
  geom_linerange(aes(ymin = Q_0.25, ymax = Q_0.75), size = 3, col = "orange") + theme_bw() +
  geom_point(size = 5, show.legend = F, shape = 23, fill = "orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_x_discrete(name = "") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "", limits = c(0,1), breaks = seq(0,1,0.1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

Figure_2 <- Figure_2A + ggtitle("(a)") + theme(plot.title = element_text(face = "bold")) +
  plot_spacer() +
  Figure_2Ba + ggtitle("(b)") + theme(plot.title = element_text(face = "bold")) +
  Figure_2Bb + ggtitle("") + theme(plot.title = element_text(face = "bold")) +
  Figure_2Bc + ggtitle("") + theme(plot.title = element_text(face = "bold")) +
  plot_spacer() +
  Figure_2C + ggtitle("(c)") + theme(plot.title = element_text(face = "bold")) +
  plot_layout(widths = c(2, 1, 8, 1, 1, 1, 10))  & 
  theme(plot.tag = element_text(face = "bold"), text = element_text(size = 20), plot.margin = unit(rep(0.1,4),"cm"))

#### Figure 3 – Simulated tradeoffs CCA vs. coral ----

# Define quantile rates
quantile(data_cca$Rate_std, probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99))
quantile(data_Niklas$conX , probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99))

# range of cover for both CCA y Coral
cca_cover       = seq(0,100,0.01)
coral_cover     = seq(100,0,-0.01)

# calcification rates for CCA and corals (kg CaCO3 m^-2 yr_1) excluding 3 outliers in CCA data
cca_rate        = quantile(data_cca$Rate_std, probs = 0.50) * 10
cca_rate_lwr    = quantile(data_cca$Rate_std, probs = 0.25) * 10
cca_rate_upr    = quantile(data_cca$Rate_std, probs = 0.75) * 10
coral_rate      = quantile(data_Niklas$conX , probs = 0.50) * 10
coral_rate_lwr  = quantile(data_Niklas$conX , probs = 0.25) * 10
coral_rate_upr  = quantile(data_Niklas$conX , probs = 0.75) * 10

# coral rugosity = 1.5, 3 assuming a flat reef rugosity = 1 for CCA
r_branching     = 2.97
r_massive       = 1.54

## Calculations for both CCA and Corals
lwr             = cca_rate_lwr * ((100 - coral_cover)/100) + coral_rate_upr * (coral_cover/100)
median          = cca_rate * ((100 - coral_cover)/100) + coral_rate * (coral_cover/100)
upr             = cca_rate_upr * ((100 - coral_cover)/100) + coral_rate_lwr * (coral_cover/100)
lwr_ccoral      = coral_rate_upr * (coral_cover/100)
median_ccoral   = coral_rate * (coral_cover/100)
upr_ccoral      = coral_rate_lwr * (coral_cover/100)
lwr_ccca        = cca_rate_upr * (cca_cover/100)
median_ccca     = cca_rate * (cca_cover/100)
upr_ccca        = cca_rate_lwr * (cca_cover/100)

lwr_b           = cca_rate_lwr * ((100 - coral_cover)/100) + coral_rate_upr * r_branching * (coral_cover/100) 
median_b        = cca_rate * ((100 - coral_cover)/100) + coral_rate * r_branching * (coral_cover/100)
upr_b           = cca_rate_upr * ((100 - coral_cover)/100) + coral_rate_lwr * r_branching * (coral_cover/100)
lwr_ccoral_b    = coral_rate_upr * r_branching * (coral_cover/100)
median_ccoral_b = coral_rate * r_branching * (coral_cover/100)
upr_ccoral_b    = coral_rate_lwr * r_branching * (coral_cover/100)

lwr_m           = cca_rate_lwr * ((100 - coral_cover)/100) + coral_rate_upr * r_massive * (coral_cover/100) 
median_m        = cca_rate * ((100 - coral_cover)/100) + coral_rate * r_massive * (coral_cover/100)
upr_m           = cca_rate_upr * ((100 - coral_cover)/100) + coral_rate_lwr * r_massive * (coral_cover/100)

lwr_ccoral_m    = coral_rate_upr * r_massive * (coral_cover/100)
median_ccoral_m = coral_rate * r_massive * (coral_cover/100)
upr_ccoral_m    = coral_rate_lwr * r_massive * (coral_cover/100)

# Generate overal data for Figure 3 
data            <- as.data.frame(cbind(cca_cover,coral_cover)) %>%
  mutate(G_percent_coral_lwr = (lwr_ccoral / lwr)*100) %>%
  mutate(G_percent_coral = (median_ccoral / median)*100) %>%
  mutate(G_percent_coral_upr = (upr_ccoral / upr)*100) %>%
  mutate(G_percent_cca_lwr = (lwr_ccca / upr)*100) %>%
  mutate(G_percent_cca = (median_ccca / median)*100) %>%
  mutate(G_percent_cca_upr = (upr_ccca / lwr)*100) %>%
  mutate(G_percent_branching_lwr = (lwr_ccoral_b / lwr_b)*100) %>%
  mutate(G_percent_branching = (median_ccoral_b / median_b)*100) %>%
  mutate(G_percent_branching_upr = (upr_ccoral_b / upr_b)*100) %>%
  mutate(G_percent_massive_lwr = (lwr_ccoral_m / lwr_m)*100) %>%
  mutate(G_percent_massive = (median_ccoral_m / median_m)*100) %>%
  mutate(G_percent_massive_upr = (upr_ccoral_m / upr_m)*100) %>%
  mutate(Gb_percent_cca_lwr = (lwr_ccca / upr_b)*100) %>%
  mutate(Gb_percent_cca = (median_ccca / median_b)*100) %>%
  mutate(Gb_percent_cca_upr = (upr_ccca / lwr_b)*100) %>%
  mutate(Gm_percent_cca_lwr = (lwr_ccca / upr_m)*100) %>%
  mutate(Gm_percent_cca = (median_ccca / median_m)*100) %>%
  mutate(Gm_percent_cca_upr = (upr_ccca / lwr_m)*100)

# Calcul treshold
flat_threshold      = data$coral_cover[ceiling(data$G_percent_coral)==50][1]
massive_threshold   = data$coral_cover[ceiling(data$G_percent_massive)==50][1]
branching_threshold = data$coral_cover[ceiling(data$G_percent_branching)==50][1]

## Rugosity treshold
# Quantify threshold tipping points between CCA and coral calcification

for(i in 1:400){
  median_r_lwr = cca_rate_lwr * ((100 - coral_cover)/100) + coral_rate_upr * (i/100) * (coral_cover/100)
  median_r     = cca_rate * ((100 - coral_cover)/100) + coral_rate * (i/100) * (coral_cover/100)
  median_r_upr = cca_rate_upr * ((100 - coral_cover)/100) + coral_rate_lwr * (i/100) * (coral_cover/100)
  
  median_ccoral_r_lwr = coral_rate_upr * (i/100) * (coral_cover/100)
  median_ccoral_r     = coral_rate * (i/100) * (coral_cover/100)
  median_ccoral_r_upr = coral_rate_lwr * (i/100) * (coral_cover/100)
  
  data_r <- data %>%
    mutate(G_percent_rugosity_lwr = (median_ccoral_r_lwr / median_r_lwr)*100) %>%
    mutate(G_percent_rugosity = (median_ccoral_r / median_r)*100) %>%
    mutate(G_percent_rugosity_upr = (median_ccoral_r_upr / median_r_upr)*100)
  
  r_lwr_threshold = data_r$coral_cover[ceiling(data_r$G_percent_rugosity_lwr) == 50][1]
  r_threshold     = data_r$coral_cover[ceiling(data_r$G_percent_rugosity) == 50][1]
  r_upr_threshold = data_r$coral_cover[ceiling(data_r$G_percent_rugosity_upr) == 50][1]
  
  rugosity.list[[i]] <- data.frame(r = i/100, 
                                   threshold_lwr = r_lwr_threshold,
                                   threshold     = r_threshold,
                                   threshold_upr = r_upr_threshold)
}

# Starting to plot Figure 3 – Surface Area Normalized Calcification rates

rates_plot_df <- data.frame(coral_cover     = coral_cover,
                            cca_cover       = cca_cover,
                            median_ccoral   = median_ccoral,
                            lwr_ccoral      = lwr_ccoral,
                            upr_ccoral      = upr_ccoral,
                            median_ccoral_b = median_ccoral_b,
                            lwr_ccoral_b    = lwr_ccoral_b,
                            upr_ccoral_b    = upr_ccoral_b,
                            median_ccoral_m = median_ccoral_m,
                            lwr_ccoral_m    = lwr_ccoral_m,
                            upr_ccoral_m    = upr_ccoral_m,
                            median_ccca     = median_ccca,
                            lwr_ccca        = lwr_ccca,
                            upr_ccca        = upr_ccca)

rates_plot_manual_legend <- data.frame(x = c(rep(30,4),rep(25,4)),
                                       y = c(seq(14,11),seq(14,11)),
                                       colour = c(rep(seq(1,4), 2))) %>% mutate(colour = factor(colour))

Figure_3A = ggplot(data=rates_plot_df) +
  geom_ribbon(aes(x=coral_cover, ymin=lwr_ccoral_b,ymax=upr_ccoral_b), fill = "#a1dab4",alpha=0.5) +
  geom_ribbon(aes(x=coral_cover, ymin=lwr_ccoral_m,ymax=upr_ccoral_m), fill = "#fee090",alpha=0.5) +
  geom_ribbon(aes(x=coral_cover, ymin=lwr_ccoral,ymax=upr_ccoral), fill = "#f46d43",alpha=0.5) +
  geom_smooth(aes(x=coral_cover,y=median_ccoral_b), color="#a1dab4",size=1.5) +
  geom_smooth(aes(x=coral_cover,y=median_ccoral_m), color="#fee090",size=1.5) +
  geom_smooth(aes(x=coral_cover,y=median_ccoral), color="#f46d43",size=1.5) +
  xlab("% Coral Cover") +
  ylab(expression(Coral~Carbonate~Production~(kg~CaCO[3]~m^-2~yr^-1))) +
  geom_smooth(aes(x=rev(cca_cover),y=median_ccca),color="#41b6c4",size=1.5) +
  geom_ribbon(aes(x=rev(cca_cover), ymin=lwr_ccca,ymax=upr_ccca), fill = "#41b6c4",alpha=0.5) +
  geom_textpath(data = rates_plot_df %>% filter(between(coral_cover,90,98)), 
                aes(x=coral_cover, y=lwr_ccoral_b-0.5),label = "Branching Coral", text_only = TRUE, colour = "black", size = 5, hjust = 0) + 
  geom_textpath(data = rates_plot_df %>% filter(between(coral_cover,90,98)), 
                aes(x=coral_cover, y=lwr_ccoral_m-0.5),label = "Massive Coral", text_only = TRUE, colour = "black", size = 5, hjust = 0) + 
  geom_textpath(data = rates_plot_df %>% filter(between(coral_cover,90,98)),
                aes(x=coral_cover, y=lwr_ccoral-0.5),label = "Coral", text_only = TRUE, colour = "black", size = 5, hjust = 0) + 
  geom_textpath(data = rates_plot_df %>% filter(between(cca_cover,78,95)), 
                aes(x=100-cca_cover, y=lwr_ccca-0.5),label = "CCA", text_only = TRUE, colour = "black", size = 5, hjust = 0) + 
  geom_line(data = rates_plot_manual_legend, aes(x = x, y = y, colour = colour), size = 3) +
  annotate(geom = "text", x = 22.5, y = seq(14,11), hjust = 0, size = 5, label = c(expression(paste(R[Coral],' = ',2.97)), 
                                                                                   expression(paste(R[Coral],' = ',1.54)), 
                                                                                   expression(paste(R[Coral],' = ',1,'.00')), 
                                                                                   expression(paste(R[CCA],' = ',1,'.00')))) +
  scale_y_continuous(name=expression(Coral~G~(kg~CaCO[3]~m^-2~yr^-1)),sec.axis=sec_axis(~.,name=expression(CCA~G~(kg~CaCO[3]~m^-2~yr^-1))), expand=c(0,0)) +
  scale_x_continuous(name="% Coral Cover",trans="reverse",sec.axis=sec_axis(~100-.,name="% CCA Cover"), expand=c(0,0)) +
  scale_colour_manual(values = rates_plot_manual_legend_colours) +
  coord_cartesian(ylim=c(0,15)) +
  ggtitle("(a) Surface Area Normalised Calcification Rates") +
  theme_classic() +
  theme(text = element_text(size=22),
        plot.title = element_text(size=17,face="bold"),
        axis.title.x.top = element_text(colour = "#41b6c4",size=16),
        axis.title.x = element_text(colour = "black",size=16),
        axis.title.y = element_text(colour = "black",size=16),
        axis.title.y.right = element_text(colour = "#41b6c4",size=16),
        axis.text.x.top = element_text(colour = "#41b6c4"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.y.right = element_text(colour = "#41b6c4"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(), 
        legend.position = "none")

# Starting to plot Figure 3 – Rugosity tresholds

Figure_3B = data.table::rbindlist(rugosity.list) %>%
  ggplot(., aes(x = r)) + 
  geom_line(aes(y = 100-threshold),size=1.5) +
  geom_ribbon(aes(x=r, ymin=100-threshold_lwr,ymax=100-threshold_upr), fill = "black",alpha=0.1) +
  geom_vline(xintercept=1,size=1.5,color="#f46d43") +
  geom_vline(xintercept=1.54,size=1.5,color="#fee090") +
  geom_vline(xintercept=2.97,size=1.5,color="#a1dab4") + 
  geom_textpath(x=2.97+0.1, y=5, label = "Branching Coral : CCA", text_only = TRUE, colour = "black", size = 4.5, hjust = 1, angle = -90) +
  geom_textpath(x=1.54+0.1, y=5, label = "Massive Coral : CCA", text_only = TRUE, colour = "black", size = 4.5, hjust = 1, angle = -90) +
  geom_textpath(x=1+0.1, y=5, label = "Coral : CCA", text_only = TRUE, colour = "black", size = 4.5, hjust = 1, angle = -90) +
  xlab(expression(paste(R[Coral]:R[CCA]))) +
  ylab(expression(paste('% CCA Cover at ',G[CCA],'>',G[Coral]))) +
  coord_cartesian(xlim=c(0,4),ylim=c(0,100)) +
  scale_y_continuous(expand=c(0,0), sec.axis = sec_axis(~100-., name = expression(paste('% Coral Cover at ', G[CCA],'>',G[Coral])))) +
  scale_x_continuous(expand=c(0,0)) +
  ggtitle(expression(bold(paste('(b) Effect of ',R[Coral]:R[CCA],' on Threshold for ', G[CCA],'>',G[Coral])))) +
  theme_classic() +
  theme(text = element_text(size=20),
        plot.title = element_text(size=17),
        axis.title.x = element_text(colour = "black",size=16),
        axis.title.y = element_text(colour = "#41b6c4",size=16),
        axis.title.y.right = element_text(colour = "black",size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "#41b6c4"),
        axis.text.y.right = element_text(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank())

# Starting to plot Figure 3 – Coral roughness = 2.97

Figure_3C = ggplot(data=data) +
  geom_line(aes(x=coral_cover,y=G_percent_branching), color="#a1dab4",size=1.5) +
  geom_ribbon(aes(x=coral_cover, ymin=G_percent_branching_lwr,ymax=G_percent_branching_upr), fill = "#a1dab4",alpha=0.5) +
  xlab("% Coral Cover") +
  ylab("% Contribution by Coral") +
  geom_line(aes(x=rev(cca_cover),y=Gb_percent_cca),color="#41b6c4",size=1.5) +
  geom_ribbon(aes(x=rev(cca_cover), ymin=Gb_percent_cca_lwr,ymax=Gb_percent_cca_upr), fill = "#41b6c4",alpha=0.5) +
  scale_y_continuous(name="% Contribution by Coral",sec.axis=sec_axis(~.,name="% Contribution by CCA"), expand=c(0,0)) +
  scale_x_continuous(name="% Coral Cover",trans="reverse",sec.axis=sec_axis(~100-.,name="% CCA Cover"), expand=c(0,0)) +
  geom_image(data=branching_svg,mapping=aes(x,y,image=image),size=0.3) +
  coord_cartesian(ylim=c(0,100)) +
  ggtitle(expression(bold(paste('(c) ',R[Coral]:R[CCA],' = ','2.97')))) +
  geom_vline(xintercept=branching_threshold,alpha=0.5,size=1.5) +
  theme_classic() +
  theme(text = element_text(size=22),
        plot.title = element_text(size=17),
        axis.title.x.top = element_text(colour = "#41b6c4",size=16),
        axis.title.x = element_text(colour = "black",size=16),
        axis.title.y = element_text(colour = "black",size=16),
        axis.title.y.right = element_text(colour = "#41b6c4",size=16),
        axis.text.x.top = element_text(colour = "#41b6c4"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.y.right = element_text(colour = "#41b6c4"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank()) +
  geom_textpath(data = data %>% filter(between(coral_cover,50,70)), 
                aes(y = G_percent_branching-5, x = coral_cover), label = "Branching Coral", text_only = TRUE, size = 5, hjust = 0.5, colour = "black") + 
  geom_textpath(data = data %>% filter(between(coral_cover,60,70)), 
                aes(y = Gb_percent_cca+5, x = coral_cover), label = "CCA", text_only = TRUE, size = 5, hjust = 0.5, colour = "black")

# Starting to plot Figure 3 – Coral roughness = 1.54

Figure_3D = ggplot(data=data) +
  geom_line(aes(x=coral_cover,y=G_percent_massive), color="#fee090",size=1.5) +
  geom_ribbon(aes(x=coral_cover, ymin=G_percent_massive_lwr,ymax=G_percent_massive_upr), fill = "#fee090",alpha=0.5) +
  xlab("% Coral Cover") +
  ylab("% Contribution by Coral") +
  geom_line(aes(x=rev(cca_cover),y=Gm_percent_cca),color="#41b6c4",size=1.5) +
  geom_ribbon(aes(x=rev(cca_cover), ymin=Gm_percent_cca_lwr,ymax=Gm_percent_cca_upr), fill = "#41b6c4",alpha=0.5) +
  scale_y_continuous(name="% Contribution by Coral",sec.axis=sec_axis(~.,name="% Contribution by CCA"), expand=c(0,0)) +
  scale_x_continuous(name="% Coral Cover",trans="reverse",sec.axis=sec_axis(~100-.,name="% CCA Cover"), expand=c(0,0)) +
  geom_image(data=massive_svg,mapping=aes(x,y,image=image),size=0.3) +
  coord_cartesian(ylim=c(0,100)) +
  ggtitle(expression(bold(paste('(d) ',R[Coral]:R[CCA],' = ','1.54')))) +
  geom_vline(xintercept=massive_threshold,alpha=0.5,size=1.5) +
  theme_classic() +
  theme(text = element_text(size=22),
        plot.title = element_text(size=17),
        axis.title.x.top = element_text(colour = "#41b6c4",size=16),
        axis.title.x = element_text(colour = "black",size=16),
        axis.title.y = element_text(colour = "black",size=16),
        axis.title.y.right = element_text(colour = "#41b6c4",size=16),
        axis.text.x.top = element_text(colour = "#41b6c4"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.y.right = element_text(colour = "#41b6c4"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank()) +
  geom_textpath(data = data %>% filter(between(coral_cover,60,70)), 
                aes(y = G_percent_massive-5, x = coral_cover), label = "Massive Coral", text_only = TRUE, size = 5, hjust = 0.5, colour = "black") + 
  geom_textpath(data = data %>% filter(between(coral_cover,60,70)), 
                aes(y = Gm_percent_cca+5, x = coral_cover), label = "CCA", text_only = TRUE, size = 5, hjust = 0.5, colour = "black")

# Starting to plot Figure 3 – Coral roughness = 1.00

Figure_3E = ggplot(data=data) +
  geom_line(aes(x=coral_cover,y=G_percent_coral), color="#f46d43",size=1.5) +
  geom_ribbon(aes(x=coral_cover, ymin=G_percent_coral_lwr,ymax=G_percent_coral_upr), fill = "#f46d43",alpha=0.5) +
  xlab("% Coral Cover") +
  ylab("% Contribution by Coral") +
  geom_line(aes(x=rev(cca_cover),y=G_percent_cca),color="#41b6c4",size=1.5) +
  geom_ribbon(aes(x=rev(cca_cover), ymin=G_percent_cca_lwr,ymax=G_percent_cca_upr), fill = "#41b6c4",alpha=0.5) +
  scale_y_continuous(name="% Contribution by Coral",sec.axis=sec_axis(~.,name="% Contribution by CCA"), expand=c(0,0)) +
  scale_x_continuous(name="% Coral Cover",trans="reverse",sec.axis=sec_axis(~100-.,name="% CCA Cover"), expand=c(0,0)) +
  geom_image(data=coral_cca_svg,mapping=aes(x,y,image=image),size=0.3) +
  coord_cartesian(ylim=c(0,100)) +
  ggtitle(expression(bold(paste('(e) ',R[Coral]:R[CCA],' = ','1.00')))) +
  geom_vline(xintercept=flat_threshold,alpha=0.5,size=1.5) +
  theme_classic() +
  theme(text = element_text(size=22),
        plot.title = element_text(size=17),
        axis.title.x.top = element_text(colour = "#41b6c4",size=16),
        axis.title.x = element_text(colour = "black",size=16),
        axis.title.y = element_text(colour = "black",size=16),
        axis.title.y.right = element_text(colour = "#41b6c4",size=16),
        axis.text.x.top = element_text(colour = "#41b6c4"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.y.right = element_text(colour = "#41b6c4"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
  geom_textpath(data = data %>% filter(between(coral_cover,60,70)), 
                aes(y = G_percent_coral-5, x = coral_cover), label = "Coral", text_only = TRUE, size = 5, hjust = 0.5, colour = "black") + 
  geom_textpath(data = data %>% filter(between(coral_cover,60,70)), 
                aes(y = G_percent_cca+5, x = coral_cover), label = "CCA", text_only = TRUE, size = 5, hjust = 0.5, colour = "black")

# Combine Figure 3
Figure_3 = (Figure_3A|Figure_3B) / (Figure_3C | Figure_3D | Figure_3E ) & 
  theme(plot.margin = unit(c(0.2,0.3,0.2,0.2), units = "in"),
        axis.text   = element_text(size = 14))

#### Figure 4 – Case Study Moorea ----

# Extract CCA cover information from LTER
LTER_cca <- LTER %>% dplyr::filter(., Taxonomy_Substrate_Functional_Group == "Crustose Corallines")

# Extract coral carbonate calcification using buoyant weight
corals_all = data_Niklas %>% group_by(Genus) %>% summarise(mean = mean(conX), sd = sd(conX))
corals_avg = data.frame(Genus = "All Corals", data_Niklas %>% summarise(mean = mean(conX), sd = sd(conX)))
corals_all = rbind(corals_all, corals_avg)

# Structural complexity from Carlot et al. (2023)
SC = data.frame(Year = c(2005, 2008:2016), SC = c(3.07, 2.73, 2.46, 1.75, 1.58, 2.21, 1.99, 2.52, 2.30, 3.87))

## Datasets management
# CCA
LTER_cca <- LTER %>% dplyr::filter(., Taxonomy_Substrate_Functional_Group == "Crustose Corallines") %>% 
  dplyr::filter(., Year %in% SC$Year)
LTER_Summary_cca <- LTER_cca %>% group_by(Year, Site, Habitat, Transect, Quadrat) %>% 
  summarise(CCA_Cover = mean(Percent_Cover)) %>% group_by(Year, Site, Habitat, Transect) %>% 
  summarise(CCA_Cover = mean(CCA_Cover)) %>% group_by(Year, Site, Habitat) %>% 
  summarise(sd_CCA_Cover = sd(CCA_Cover), CCA_Cover = mean(CCA_Cover)) %>% group_by(Year, Site) %>% 
  summarise(CCA_Cover = mean(CCA_Cover), sd_CCA_Cover = mean(sd_CCA_Cover)) %>% 
  mutate(., CR = CCA_Cover/100 * cca_rate)
LTER_CCA_avg = LTER_Summary_cca %>% group_by(Year) %>% 
  summarise(Cover = mean(CCA_Cover), sd_cover = sd(CCA_Cover), sd_CR = sd(CR), CR = mean(CR)) %>% 
  mutate(., Cover = as.numeric(Cover), sd_cover = as.numeric(sd_cover), CR = as.numeric(CR),
         sd_CR = as.numeric(sd_CR), Year = as.numeric(Year))

# Coral
LTER_coral <- LTER %>% dplyr::filter(., Taxonomy_Substrate_Functional_Group == "Coral") %>% 
  dplyr::filter(., Year %in% SC$Year)
LTER_Summary_coral <- LTER_coral %>% group_by(Year, Site, Habitat, Transect, Quadrat) %>% 
  summarise(CCA_Cover = mean(Percent_Cover)) %>% group_by(Year, Site, Habitat, Transect) %>% 
  summarise(CCA_Cover = mean(CCA_Cover)) %>% group_by(Year, Site, Habitat) %>% 
  summarise(sd_CCA_Cover = sd(CCA_Cover), CCA_Cover = mean(CCA_Cover)) %>% group_by(Year, Site) %>% 
  summarise(CCA_Cover = mean(CCA_Cover), sd_CCA_Cover = mean(sd_CCA_Cover)) %>% 
  mutate(., CR = CCA_Cover/100 * coral_rate) %>% mutate(., Year = as.numeric(Year)) %>% 
  inner_join(., SC, by = "Year") %>% mutate(., CR_SC = CR * SC)
LTER_Coral_avg = LTER_Summary_coral %>% group_by(Year) %>% 
  summarise(sd_cover = sd(CCA_Cover), Cover = mean(CCA_Cover), sd_CR = sd(CR), CR = mean(CR), CR = mean(CR)) %>% 
  mutate(., Cover = as.numeric(Cover), sd_cover = as.numeric(sd_cover), CR = as.numeric(CR),
         sd_CR = as.numeric(sd_CR), Year = as.numeric(Year), CR_SC = CR * SC$SC, sd_CR_SC = sd_CR * SC$SC)

# Contribution definition
Contribution = data.frame(Year = LTER_Summary_cca$Year, 
                          Site = LTER_Summary_cca$Site,
                          Contribution = (LTER_Summary_cca$CR / (LTER_Summary_cca$CR + LTER_Summary_coral$CR_SC)) * 100)
Contribution_avg = Contribution %>% group_by(Year) %>% 
  summarise(Contribution_avg = mean(Contribution), sd = sd(Contribution)) %>% 
  mutate(., upr = Contribution_avg+sd, lwr = Contribution_avg-sd) %>% 
  mutate(., Year = as.numeric(Year), upr = as.numeric(upr), lwr = as.numeric(lwr))

## Starting to plot Figure 4 – Case study
Figure_4A <- ggplot(LTER_CCA_avg, aes(x = Year - 0.2, y = Cover)) + 
  geom_line(col = "pink", linetype = "dashed", size = .5) +
  geom_linerange(aes(ymin = Cover - sd_cover, ymax = Cover + sd_cover), col = "violetred") + theme_bw() +
  geom_line(data = LTER_Coral_avg, aes(x = Year + 0.2, y = Cover), col = "gold", linetype = "dashed", size = .5) +
  geom_linerange(data = LTER_Coral_avg, aes(x = Year + 0.2, ymin = Cover - sd_cover, ymax = Cover + sd_cover), col = "orange") +
  geom_point(size = 4, show.legend = F, shape = 21, fill = "violetred") +
  geom_point(data = LTER_Coral_avg, aes(x = Year + 0.2, y = Cover), size = 4, show.legend = F, shape = 23, fill = "orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "", breaks = seq(2005, 2020,5)) +
  scale_y_continuous(name = "Cover (%)", limits = c(0,50), breaks = seq(0,50,10)) 

Figure_4B <- ggplot(LTER_CCA_avg, aes(x = Year - 0.2, y = CR)) + 
  geom_line(col = "pink", linetype = "dashed", size = .5) +
  geom_linerange(aes(ymin = CR - sd_CR, ymax = CR + sd_CR), col = "violetred") + theme_bw() +
  geom_line(data = LTER_Coral_avg, aes(x = Year + 0.2, y = CR_SC), col = "gold", linetype = "dashed", size = .5) +
  geom_linerange(data = LTER_Coral_avg, aes(x = Year + 0.2, ymin = CR_SC - sd_CR_SC, ymax = CR_SC + sd_CR_SC), col = "orange") +
  geom_point(size = 4, show.legend = F, shape = 21, fill = "violetred") +
  geom_point(data = LTER_Coral_avg, aes(x = Year + 0.2, y = CR_SC), size = 4, show.legend = F, shape = 23, fill = "orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "", breaks = seq(2005, 2020,5)) + 
  scale_y_continuous(name = expression(atop(CaCO[3]~"production", paste("(kg."*m^-2*".yr"^-1*")"))), 
                     limits = c(0,6), breaks = seq(0,6,1)) 

Figure_4C <- ggplot(Contribution_avg, aes(x = as.numeric(Year), y = as.numeric(Contribution_avg))) + 
  geom_line(col = "pink", linetype = "dashed", size = .5) +
  geom_linerange(aes(ymin = as.numeric(Contribution_avg) - as.numeric(sd), 
                     ymax = as.numeric(Contribution_avg) + as.numeric(sd)), col = "violetred") + theme_bw() +
  geom_point(size = 4, show.legend = F, shape = 21, fill = "violetred") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "", breaks = seq(2005,2020,5)) + 
  scale_y_continuous(name = expression(atop("CCA contribution to", paste("the"~CaCO[3]~"budget (%)"))), 
                     limits = c(0,100), breaks = seq(0,100,10)) 

Figure_4 <- Figure_4A + ggtitle("(a)") + theme(plot.title = element_text(face = "bold")) +
  Figure_4B + ggtitle("(b)") + theme(plot.title = element_text(face = "bold")) +
  Figure_4C + ggtitle("(c)") + theme(plot.title = element_text(face = "bold"))  & 
  theme(plot.tag = element_text(face = "bold"), text = element_text(size = 20), plot.margin = unit(rep(0.5,4),"cm")) &
  scale_shape_manual(name = "", values = c(21,23),                 limits = c("CCA", "Corals")) &
  scale_color_manual(name = "", values = c("violetred", "orange"), limits = c("CCA", "Corals")) 

#### Figure 5 – Systematic review vizualisation ----

{GR_viz[[1]] <- Data_viz[[1]] %>% dplyr::filter(., `Paper name` %notin% c("Graba-Landry et al. 2018", "Johnson et al 2019"),
                                               Genus != "Halimeda") %>%
  dplyr::filter(., Genus %notin% data_articulate$Genus) %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std - std_error, ymax = Rate_std + std_error, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = expression("Calcification rate (µmol."*g^-1*".h"^-1*")"), limits = c(0,30))}
{GR_viz[[2]] <- Data_viz[[1]] %>% dplyr::filter(., `Paper name` %notin% c("Graba-Landry et al. 2018", "Johnson et al 2019"),
                                               Genus != "Halimeda") %>%
  dplyr::filter(., Genus %in% data_articulate$Genus) %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std - std_error, ymax = Rate_std + std_error, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + scale_y_continuous(name = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())}
{GR_viz[[3]] <- Data_viz[[1]] %>% dplyr::filter(., `Paper name` %notin% c("Graba-Landry et al. 2018", "Johnson et al 2019"),
                                               Genus == "Halimeda") %>%
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std - std_error, ymax = Rate_std + std_error, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + scale_y_continuous(name = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())}
{GR_viz[[4]] <- Data_viz[[2]] %>% dplyr::filter(., Genus %notin% c("Halimeda")) %>% 
  dplyr::filter(., Genus %notin% data_articulate$Genus) %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std - std_error, ymax = Rate_std + std_error, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + scale_y_continuous(name = expression("Linear extension (mm."*yr^-1*")")) }
{GR_viz[[5]] <- Data_viz[[2]] %>% dplyr::filter(., Genus %notin% c("Halimeda")) %>% 
  dplyr::filter(., Genus %in% data_articulate$Genus) %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std - std_error, ymax = Rate_std + std_error, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + scale_y_continuous(name = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())}
{GR_viz[[6]] <- Data_viz[[3]] %>% dplyr::filter(., Genus %notin% c("Halimeda")) %>% 
  dplyr::filter(., Genus %notin% data_articulate$Genus) %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std * 0.36525, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std * 0.36525 - std_error * 0.36525, 
                     ymax = Rate_std * 0.36525 + std_error * 0.36525, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = expression("Calcification rate (g."*cm^-2*".yr"^-1*")"), 
                     limits = c(0,1.2), breaks = seq(0,1.2,.3))}
{GR_viz[[7]] <- Data_viz[[3]] %>% dplyr::filter(., Genus %notin% c("Halimeda")) %>% 
  dplyr::filter(., Genus %in% data_articulate$Genus) %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std * 0.36525, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std * 0.36525 - std_error * 0.36525, 
                     ymax = Rate_std * 0.36525 + std_error * 0.36525, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "", limits = c(0,1.2), breaks = seq(0,1.2,.3)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())}
{GR_viz[[8]] <- Data_viz[[3]] %>% dplyr::filter(., Genus == "Halimeda") %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std * 0.36525, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std * 0.36525 - std_error * 0.36525, 
                     ymax = Rate_std * 0.36525 + std_error * 0.36525, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "", limits = c(0,1.2), breaks = seq(0,1.2,.3)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())}

## Generating a mixed bayesian model using BRMS to unravel what's the best technique to use?
# No enough values for CT scan – We cannot use those observations for the model
Data_viz_1 = Data_viz[[1]] %>% dplyr::filter(., Method_family != "X-ray CT Scan", 
                                             `Paper name` %notin% c("Graba-Landry et al. 2018", "Johnson et al 2019"),
                                             Genus %notin% c("Halimeda"))
Data_viz_2 = Data_viz[[2]] %>% dplyr::filter(., Genus %notin% c("Halimeda"))
Data_viz_3 = Data_viz[[3]] %>% dplyr::filter(., Genus %notin% c("Halimeda"))

# To gain time, models have been saved into models folder
# fit_biom      <- brms::brm(Rate_std ~ Climate + Genus + (1 | Method_family), data = Data_viz_1, family = gaussian(), 
#                            warmup = 3000, iter = 5000, control = list(max_treedepth = 10, adapt_delta = 0.99), core = 4)
# fit_extension <- brms::brm(Rate_std ~ Climate + Genus + (1 | Method_family), data = Data_viz_2, family = gaussian(), 
#                            warmup = 3000, iter = 5000, control = list(max_treedepth = 10, adapt_delta = 0.99), core = 4)
# fit_surface   <- brms::brm(Rate_std ~ Climate + Genus + (1 | Method_family), data = Data_viz_3, family = gaussian(), 
#                            warmup = 3000, iter = 5000, control = list(max_treedepth = 11, adapt_delta = 0.99), core = 4)}

# Effect size
summary_extension = fit_extension %>% spread_draws(r_Method_family[condition,]) %>% summarise_draws()
summary_surface   = fit_surface %>% spread_draws(r_Method_family[condition,]) %>% summarise_draws()
summary_biomass   = fit_biom %>% spread_draws(r_Method_family[condition,]) %>% summarise_draws()

# Summarizing a unique dataset
summary           = rbind(summary_biomass, summary_extension, summary_surface)
summary$Methods   = chartr(".", " ", summary$condition)
for (i in 1:8) { vector_method[[i]] <- rnorm(100000, summary$mean[i], summary$sd[i]) }
data_methods = data.frame(Methods = rep(c("BW or RGR", "Isotopes", "TA anomaly", "BW or RGR", "Staining", 
                                          "BW or RGR", "TA anomaly", "X-ray CT Scan"), each = 100000),
                          Random_effect = abind::abind(vector_method), 
                          dataset = c(rep("Biomass", 300000), rep("Extension", 200000), rep("Surface", 300000))) 

# Re-arange the data in 3 object regarding the metric used
data_methods = complete(data_methods, Methods, dataset, fill = list(Random_effect = NA)) %>% group_split(dataset)
summary_list[[1]] = summary[1:3,] 
summary_list[[2]] = summary[4:5,] 
summary_list[[3]] = summary[6:8,]

# Re-arrange the colors for each dataset
color_list[[1]] = colors[c(1:2,4)] 
color_list[[2]] = colors[c(1,3)] 
color_list[[3]] = colors[c(1,4:5)]

# Starting to plot Figure 5 – First part: systematic review viz
Figure_5A <-  GR_viz[[1]] + ggtitle("(a)") + theme(plot.title = element_text(face = "bold")) +
              GR_viz[[2]] + ggtitle("   ") + theme(plot.title = element_text(face = "bold")) +
              GR_viz[[3]] + ggtitle("   ") + theme(plot.title = element_text(face = "bold")) +
              plot_spacer() +
              GR_viz[[6]] + ggtitle("(c)") + theme(plot.title = element_text(face = "bold")) +
              GR_viz[[7]] + ggtitle("   ") + theme(plot.title = element_text(face = "bold")) +
              GR_viz[[8]] + ggtitle("   ") + theme(plot.title = element_text(face = "bold")) +
              plot_layout(guides = "collect", nrow = 1, widths = c(9,4,1,2,9,1,1)) &
              scale_shape_manual(name = "Methods", limits = c("Isotopes", "BW or RGR", "TA anomaly", "X-ray CT Scan"), 
                                 values = c(22, 21, 23, 24)) &
              scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) &
              theme(plot.margin = unit(rep(0.1,4),"cm"))

# Re-arrange dataset for Effect size
summary_list[[1]] = rbind(summary_list[[1]],
                          data.frame(condition = "X-ray.CT.Scan", variable = "r_Method_family",
                                     mean = NA, median = NA, sd = NA, mad = NA, q5 = NA,
                                     q95 = NA, rhat = NA, ess_bulk = NA, ess_tail = NA,
                                     Methods = "X-ray CT Scan"))
summary_list[[3]] = rbind(summary_list[[3]],
                          data.frame(condition = "Isotopes", variable = "r_Method_family",
                                     mean = NA, median = NA, sd = NA, mad = NA, q5 = NA,
                                     q95 = NA, rhat = NA, ess_bulk = NA, ess_tail = NA,
                                     Methods = "Isotopes"))

# Second part of the plot – Effect size
Plot_5Ba <- data_methods[[1]] %>% dplyr::filter(Methods != "Staining") %>% drop_na() %>% 
  rbind(data.frame(Methods = "X-ray CT Scan", dataset = "Biomass", Random_effect = 0)) %>% 
  ggplot(aes(y = Methods, x = Random_effect, fill = Methods, color = Methods, shape = Methods)) + 
  theme_bw() +
  geom_segment(aes(x = 0, xend = 0, y = -Inf, yend = Inf), size = 0.75, linetype = 1, color = "black") +
  geom_segment(data = summary_list[[1]], aes(x = mean-sd, xend = mean+sd, 
                                             y = Methods , yend = Methods, color = Methods), size = 2) +
  geom_point(data = summary_list[[1]], aes(x = mean, y = Methods, fill = Methods), size = 4, color = "black") +
  scale_x_continuous(name = "Effect size", limits = c(-5,5)) +
  scale_y_discrete(name = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#ff595e", "#ffca3a", "#1982c4", "#6a4c93")) + 
  scale_color_manual(values = c("#ff595e", "#ffca3a", "#1982c4", "#6a4c93")) + 
  scale_shape_manual(values=c(21, 22, 23, 24)) +
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), strip.text.x = element_text(size = 8)) + #coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Plot_5Bb <- data_methods[[3]] %>% dplyr::filter(Methods != "Staining") %>% drop_na() %>% 
  rbind(data.frame(Methods = "Isotopes", dataset = "Biomass", Random_effect = 0)) %>% 
  ggplot(aes(y = Methods, x = Random_effect, fill = Methods, color = Methods, shape = Methods)) + 
  theme_bw() +
  geom_segment(aes(x = 0, xend = 0, y = -Inf, yend = Inf), size = 1, linetype = 1, color = "black") +
  geom_segment(data = summary_list[[3]], aes(x = mean-sd, xend = mean+sd, 
                                             y = Methods , yend = Methods, color = Methods), size = 2) +
  geom_point(data = summary_list[[3]], aes(x = mean, y = Methods, fill = Methods), size = 4, color = "black") +
  scale_x_continuous(name = "Effect size", limits = c(-5,5)) +
  scale_y_discrete(name = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#ff595e", "#ffca3a", "#1982c4", "#6a4c93")) + 
  scale_color_manual(values = c("#ff595e", "#ffca3a", "#1982c4", "#6a4c93")) + 
  scale_shape_manual(values=c(21, 22, 23, 24)) +
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), strip.text.x = element_text(size = 8)) + #coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

Figure_5B <- (Plot_5Ba + ggtitle("(b)") + theme(plot.title = element_text(face = "bold", size = 15)) +
              Plot_5Bb + ggtitle("(d)") + theme(plot.title = element_text(face = "bold", size = 15))) 

#### Export Figures ----

ggsave(Figure_2, filename  = "Figures/Raw/Figure_2.png" , device = "png", height = 15, width = 35.0, units = "cm", dpi = 300)
ggsave(Figure_3, filename  = "Figures/Raw/Figure_3.png" , device = "png", height = 10, width = 14.5, units = "in", dpi = 300)
ggsave(Figure_4, filename  = "Figures/Raw/Figure_4.png" , device = "png", height = 13, width = 35.0, units = "cm", dpi = 300)
ggsave(Figure_5A, filename = "Figures/Raw/Figure_5A.png", device = "png", height = 10, width = 25.0, units = "cm", dpi = 300)
ggsave(Figure_5B, filename = "Figures/Raw/Figure_5B.png", device = "png", height = 06, width = 25.0, units = "cm", dpi = 300)