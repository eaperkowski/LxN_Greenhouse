###########################################################
## Load libraries
###########################################################
library(tidyverse)
library(ggpubr)

###########################################################
## Load libraries, create species labels, and colorblind
## friendly palette
###########################################################
source("https://raw.githubusercontent.com/eaperkowski/LxN_Greenhouse/main/R_files/create_ndemand_metrics.R")

species.label <- c("G. max", "G. hirsutum")
names(species.label) <- c("Soybean", "Cotton")

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

###########################################################
## Remove outliers per analyses (Bonferroni p<0.05)
###########################################################
df$n.cost[c(44, 87, 275, 279, 321, 337)] <- NA
df$n.acq[c(20, 290)] <- NA
df$root.carbon.mass[c(44, 87)] <- NA
df$nod.wt[c(10, 203)] <- NA

###########################################################
## Add pubtheme
###########################################################
pubtheme <- theme_bw() +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, face = "italic"),
        panel.border = element_rect(size = 3, fill = NA),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
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
  geom_segment(aes(x = 0, xend = 630, y = 4.92, yend = -4.89e-03*630 + 4.92), color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 3.55, yend = -3.31e-03*630 + 3.55), color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 2.86, yend = -2.53e-03*630 + 2.86), color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 2.06, yend = -1.68e-03*630 + 2.06), color = "#CC79A7", size = 3) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = expression(bold("Carbon cost to acquire nitrogen (gC gN"^-1~")")),
       color = "Shade cover (%)") +
  pubtheme +
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
  geom_segment(aes(x = 0, xend = 630, y = 3.50, yend = -5.01e-04*630 + 3.50), 
               color = "#E69F00", size = 3, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 630, y = 2.60, yend = -9.53e-04*630 + 2.60), color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 2.14, yend = -1.07e-03*630 + 2.14), color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 1.59, yend = -1.08e-03*630 + 1.59), color = "#CC79A7", size = 3) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5)) +
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

ggsave("/Users/eaperkowski/fig1.ncost.facet.png",
       fig1,
       height = 7,
       width = 16,
       dpi = "print")

###########################################################
## Whole-plant nitrogen mass
###########################################################
nacq.nppm.cotton <- ggplot(data = subset(df, spp == "Cotton"), 
                           aes(x = n.ppm,
                               y = n.acq,
                               color = factor(shade.cover))) +
  geom_jitter(size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 0.04, yend = 1.10e-04*630 + 0.04), color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.033, yend = 6.45e-05*630 + 0.033), color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.029, yend = 4.35e-05*630 + 0.029), color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.024, yend = 2.15e-05*630 + 0.024), color = "#CC79A7", size = 3) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 0.16), breaks = seq(0, 0.16, 0.04)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = expression(bold("Whole plant nitrogen biomass (gN)")),
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
  geom_segment(aes(x = 0, xend = 630, y = 0.055, yend = 8.91e-05*630 + 0.055), 
               color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.045, yend = 6.88e-05*630 + 0.045), 
               color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.039, yend = 5.79e-05*630 + 0.039), 
               color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.032, yend = 4.46e-05*630 + 0.032), 
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

ggsave("/Users/eaperkowski/fig2.nacq.facet.png",
       fig2,
       height = 7,
       width = 16,
       dpi = "print")

###########################################################
## Root carbon mass
###########################################################
rootcarbon.nppm.cotton <- ggplot(data = subset(df, spp == "Cotton"), 
                                 aes(x = n.ppm, 
                                     y = root.carbon.mass,
                                     color = factor(shade.cover))) +
  geom_jitter(size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 0.186, yend = 1.07e-04*630 + 0.186), 
               color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.125, yend = 4.82e-05*630 + 0.125), 
               color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.091, yend = 1.94e-05*630 + 0.091), 
               color = "#009E73", size = 3, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 630, y = 0.050, yend = -8.52e-06*630 + 0.050), 
               color = "#CC79A7", size = 3, linetype = "dashed") +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2)) +
  labs(x = "Nitrogen fertilization (ppm)",
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
  geom_segment(aes(x = 0, xend = 630, y = 0.196, yend = 2.46e-04*630 + 0.196), 
               color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.132, yend = 1.23e-04*630 + 0.132), 
               color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.097, yend = 6.20e-05*630 + 0.097), 
               color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.054, yend = -6.32e-07*630 + 0.054), 
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

ggsave("/Users/eaperkowski/fig3.rootcarbon.facet.png",
       fig3,
       height = 7,
       width = 16,
       dpi = "print")

###########################################################
## Root nodule weight
###########################################################
nod.weight.ppm <- ggplot(data = subset(df, spp == "Soybean" & nod.wt > 0), 
                           aes(x = n.ppm,
                               y = nod.wt,
                               color = factor(shade.cover))) +
  geom_jitter(size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 0.092, yend = -1.36e-04*630 + 0.092), 
               color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.062, yend = -9.66e-05*630 + 0.062), 
               color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.045, yend = -7.37e-05*630 + 0.045), 
               color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.025, yend = -4.47e-05*630 + 0.025), 
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
  scale_y_continuous(limits = c(-0.005, 0.2), breaks = seq(0, 0.2, 0.05)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = expression(bold("Root nodule biomass (g)")),
       color = "Shade cover (%)") +
  pubtheme +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA)))
nod.weight.ppm

###########################################################
## Root nodule biomass : root biomass
###########################################################
df$nod.root.ratio <- df$nod.wt / df$roots.wt

nod.root.ppm <- ggplot(data = subset(df, spp == "Soybean" & nod.root.ratio > 0),
                       aes(x = n.ppm,
                           y = nod.root.ratio,
                           color = factor(shade.cover))) +
  geom_jitter(size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 0.218, yend = -3.79e-04*630 + 0.218), 
               color = "#E69F00", size = 2) +
  geom_segment(aes(x = 0, xend = 630, y = 0.212, yend = -3.72e-04*630 + 0.212), 
               color = "#56B4E9", size = 2) +
  geom_segment(aes(x = 0, xend = 630, y = 0.208, yend = -3.68e-04*630 + 0.208), 
               color = "#009E73", size = 2) +
  geom_segment(aes(x = 0, xend = 630, y = 0.203, yend = -3.61e-04*630 + 0.203), 
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
       y = expression(bold("Root nodule biomass: Root biomass (g g"^"-1"~")")),
       color = "Shade cover (%)") +
  pubtheme +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA)))
nod.root.ppm

fig4 <- ggarrange(nod.weight.ppm, nod.root.ppm, 
               ncol = 2, common.legend = TRUE, 
               legend = "right", align = "hv", 
               labels = "AUTO", 
               font.label = list(size = 18))
fig4

ggsave("/Users/eaperkowski/fig5.nodwgt.facet.png",
       fig4,
       width = 16,
       height = 7,
       dpi = "print")
