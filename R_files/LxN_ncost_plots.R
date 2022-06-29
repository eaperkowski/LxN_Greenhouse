###########################################################
## Load libraries
###########################################################
library(tidyverse)
library(ggpubr)

###########################################################
## Load dataset, create species labels, and colorblind
## friendly palette
###########################################################
source("LxN_ncost_create_metrics.R")

species.label <- c("G. max", "G. hirsutum")
names(species.label) <- c("Soybean", "Cotton")
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

###########################################################
## Remove outliers per analyses (Bonferroni p<0.05)
###########################################################
df$n.cost[c(44, 271, 275, 316, 331, 332)] <- NA
df$n.acq[c(286)] <- NA
df$root.carbon.mass[c(44, 86)] <- NA
df$nod.wt[c(10, 203)] <- NA

###########################################################
# Remove rows with missing biomass or nitrogen data
###########################################################
df <- df %>%
  filter(complete.cases(stem.wt, leaves.wt, n.stem))

###########################################################
## Add pubtheme
###########################################################
pubtheme <- theme_bw() +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 20, face = "italic"),
        panel.border = element_rect(size = 3, fill = NA),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 24, face = "bold"),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        axis.ticks.length = unit(0.25, "cm"))

###########################################################
## Carbon cost to acquire nitrogen
###########################################################
ncost.nppm.cotton <- ggplot(data = subset(df, spp == "Cotton"), 
                     aes(x = n.ppm,
                         y = n.cost,
                         color = factor(shade.cover))) +
  geom_jitter(aes(color = factor(shade.cover)), 
              size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 4.92, yend = -4.89E-03*630 + 4.92), color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 3.55, yend = -3.31E-03*630 + 3.55), color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 2.86, yend = -2.53E-03*630 + 2.86), color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 2.06, yend = -1.68E-03*630 + 2.06), color = "#CC79A7", size = 3) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = expression(bold("Carbon cost to acquire nitrogen (gC gN"^-1~")")),
       color = "Shade cover (%)") +
  pubtheme +
  theme(axis.text.y = element_text(size = 22)) +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA))) +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
ncost.nppm.cotton

ncost.nppm.soy <- ggplot(data = subset(df, spp == "Soybean"), 
                            aes(x = n.ppm,
                                y = n.cost,
                                color = factor(shade.cover))) +
  geom_jitter(aes(color = factor(shade.cover)), 
              size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 3.52, yend = -0.000858*630 + 3.52), 
               color = "#E69F00", size = 3, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 630, y = 2.71, yend = -0.001013*630 + 2.71), 
               color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 2.23, yend = -0.001068*630 + 2.23), 
               color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 1.60, yend = -0.001078*630 + 1.60), 
               color = "#CC79A7", size = 3) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = NULL,
       color = "Shade cover (%)") +
  pubtheme +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA))) +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
ncost.nppm.soy

fig1 <- ggarrange(ncost.nppm.cotton, ncost.nppm.soy,
                  ncol = 2, common.legend = TRUE,
                  legend = "right", align = "hv", 
                  labels = NULL, 
                  font.label = list(size = 18))
fig1

###########################################################
## Whole-plant nitrogen mass
###########################################################
nacq.nppm.cotton <- ggplot(data = subset(df, spp == "Cotton"), 
                           aes(x = n.ppm,
                               y = n.acq,
                               color = factor(shade.cover))) +
  geom_jitter(size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 0.04, yend = 1.10e-04*630 + 0.04), 
               color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.033, yend = 6.46e-05*630 + 0.033), 
               color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.029, yend = 4.35e-05*630 + 0.029), 
               color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.024, yend = 2.15e-05*630 + 0.024), 
               color = "#CC79A7", size = 3) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 0.16), breaks = seq(0, 0.16, 0.04)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = "Whole plant nitrogen biomass (g N)",
       color = "Shade cover (%)") +
  pubtheme +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA))) +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
nacq.nppm.cotton

nacq.nppm.soy <- ggplot(data = subset(df, spp == "Soybean"), 
                        aes(x = n.ppm,
                            y = n.acq,
                            color = factor(shade.cover))) +
  geom_jitter(size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 0.057, yend = 8.51e-05*630 + 0.057), 
               color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.048, yend = 6.80e-05*630 + 0.048), 
               color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.042, yend = 5.77e-05*630 + 0.042), 
               color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.034, yend = 4.37e-05*630 + 0.034), 
               color = "#CC79A7", size = 3) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 0.16), breaks = seq(0, 0.16, 0.04)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = NULL,
       color = "Shade cover (%)") +
  pubtheme +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA))) +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
nacq.nppm.soy

fig2 <- ggarrange(nacq.nppm.cotton, nacq.nppm.soy,
                  ncol = 2, common.legend = TRUE,
                  legend = "right", align = "hv", 
                  labels = NULL, 
                  font.label = list(size = 18))
fig2

###########################################################
## Root carbon mass
###########################################################
rootcarbon.nppm.cotton <- ggplot(data = subset(df, spp == "Cotton"), 
                                 aes(x = n.ppm, 
                                     y = root.carbon.mass,
                                     color = factor(shade.cover))) +
  geom_jitter(size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 0.186, yend = 1.05E-04*630 + 0.186), 
               color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.125, yend = 4.75E-05*630 + 0.125), 
               color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.091, yend = 1.92E-05*630 + 0.091), 
               color = "#009E73", size = 3, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 630, y = 0.050, yend = -8.26E-06*630 + 0.050), 
               color = "#CC79A7", size = 3, linetype = "dashed") +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2)) +
  labs(x = expression(bold("Nitrogen fertilization (ppm)")),
       y = expression(bold("Root carbon biomass (gC)")),
       color = "Shade cover (%)") +
  pubtheme +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA))) +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
rootcarbon.nppm.cotton

rootcarbon.nppm.soy <- ggplot(data = subset(df, spp == "Soybean"), 
                              aes(x = n.ppm, 
                                  y = root.carbon.mass,
                                  color = factor(shade.cover))) +
  geom_jitter(size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 0.196, yend = 2.49E-04*630 + 0.192), 
               color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.132, yend = 1.25E-04*630 + 0.131), 
               color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.097, yend = 6.22E-05*630 + 0.097), 
               color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.054, yend = -4.66e-07*630 + 0.055), 
               color = "#CC79A7", size = 3, linetype = "dashed") +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = NULL,
       color = "Shade cover (%)") +
  pubtheme +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA))) +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
rootcarbon.nppm.soy

fig3 <- ggarrange(rootcarbon.nppm.cotton, rootcarbon.nppm.soy,
                  ncol = 2, common.legend = TRUE,
                  legend = "right", align = "hv", 
                  labels = NULL, 
                  font.label = list(size = 18))
fig3


###########################################################
## Root nodule weight
###########################################################
nod.weight.ppm <- ggplot(data = subset(df, spp == "Soybean" & nod.wt > 0), 
                           aes(x = n.ppm,
                               y = nod.wt,
                               color = factor(shade.cover))) +
  geom_jitter(size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 0.091, yend = -1.35e-04*630 + 0.090), 
               color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.061, yend = -9.55e-05*630 + 0.060), 
               color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.045, yend = -7.27e-05*630 + 0.044), 
               color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.025, yend = -4.37e-05*630 + 0.024), 
               color = "#CC79A7", size = 3) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("0", "30",
                               "50", "80")) +
  scale_shape_discrete(labels = species.label,
                       guide = guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  #scale_y_continuous(limits = c(-0.005, 0.2), breaks = seq(0, 0.2, 0.05)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = "Root nodule biomass (g)",
       color = "Shade cover (%)") +
  pubtheme +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA)))
nod.weight.ppm

###########################################################
## Root nodule biomass : root biomass
###########################################################
nod.root.ppm <- ggplot(data = subset(df, spp == "Soybean" & nod.root.ratio > 0),
                       aes(x = n.ppm,
                           y = nod.root.ratio,
                           color = factor(shade.cover))) +
  geom_jitter(size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 0.200, yend = -3.36E-04*630 + 0.200), 
               color = "#E69F00", size = 2) +
  geom_segment(aes(x = 0, xend = 630, y = 0.212, yend = -3.42e-04*630 + 0.198), 
               color = "#56B4E9", size = 2) +
  geom_segment(aes(x = 0, xend = 630, y = 0.208, yend = -3.46e-04*630 + 0.196), 
               color = "#009E73", size = 2) +
  geom_segment(aes(x = 0, xend = 630, y = 0.203, yend = -3.52e-04*630 + 0.194), 
               color = "#CC79A7", size = 2) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("0", "30",
                               "50", "80")) +
  scale_shape_discrete(labels = species.label,
                       guide = guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(-0.03, 0.8), breaks = seq(0, 0.8, 0.2)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = expression(bold("Nodule biomass : Root biomass (g g"^"-1"~")")),
       color = "Shade cover (%)") +
  pubtheme +
  theme(axis.title.y = element_text(size = 21)) +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA)))
nod.root.ppm

fig4 <- ggarrange(nod.weight.ppm, nod.root.ppm, 
               ncol = 2, common.legend = TRUE, 
               legend = "right", align = "hv", 
               labels = "AUTO", 
               font.label = list(size = 24))
fig4
