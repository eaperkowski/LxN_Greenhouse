###########################################################
## Load libraries
###########################################################
library(ggplot2)
library(ggpubr)

###########################################################
## Load libraries, create species labels, and colorblind
## friendly palette
###########################################################
source("LxN_ncost_create_metrics.R")

species.label <- c("G. max", "G. hirsutum")
names(species.label) <- c("Soybean", "Cotton")

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

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
## Calculate BVR, RMF, and total biomass
###########################################################
df <- df %>%
  mutate(total.bio = stem.wt + roots.wt + leaves.wt,
         rmf = roots.wt / total.bio,
         bvr = total.bio / 3) %>%
  filter(complete.cases(stem.wt, leaves.wt, n.stem))

###########################################################
## Remove outliers per analyses (Bonferroni p<0.05)
###########################################################
df$rmf[c(86, 90, 316, 332)] <- NA
df$rmf[c(44, 171, 271, 275, 331)] <- NA

###########################################################
# Biomass : volume ratio plots
###########################################################        
bvr.cotton <- lmer(log(bvr) ~ shade.cover * n.ppm + (1 | block), 
                   data = subset(df, spp == "Cotton"))
bvr.soy <- lmer(bvr ~ shade.cover * n.ppm + (1 | block), 
                data = subset(df, spp == "Soybean"))

bvr.soy.regline <- data.frame(emmeans(bvr.soy, ~shade.cover, "n.ppm",
                                     at = list(n.ppm = seq(0, 630, 5),
                                               shade.cover = c(0, 30, 50, 80)),
                                     type = "response")) 
bvr.cotton.regline <- data.frame(emmeans(bvr.cotton, ~shade.cover, "n.ppm",
                                      at = list(n.ppm = seq(0, 630, 5),
                                                shade.cover = c(0, 30, 50, 80)),
                                      type = "response")) 

bvr.soy <- ggplot(data = subset(df, spp == "Soybean"), 
                  aes(x = n.ppm, y = bvr)) +
  geom_hline(yintercept = 1, color = "red", size = 2) +
  geom_jitter(aes(fill = factor(shade.cover)), 
              size = 4.5, alpha = 0.75, width = 20, shape = 21) +
  geom_smooth(data = bvr.soy.regline,
              aes(y = emmean, 
                  color = factor(shade.cover)), 
              size = 4) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30", "50", "80")) +
  scale_fill_manual(values = cbbPalette,
                     labels = c("0", "30", "50", "80")) +
  scale_x_continuous(limits = c(-20, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = NULL, fill = "Shade cover (%)") +
  pubtheme +
  guides(color = "none") +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
bvr.soy

bvr.cotton <- ggplot(data = subset(df, spp == "Cotton"), 
                  aes(x = n.ppm, y = bvr)) +
  geom_hline(yintercept = 1, color = "red", size = 2) +
  geom_jitter(aes(fill = factor(shade.cover)), 
              size = 4.5, alpha = 0.75, width = 20, shape = 21) +
  geom_smooth(data = bvr.cotton.regline,
              aes(y = response, 
                  color = factor(shade.cover)), 
              size = 4) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30", "50", "80")) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("0", "30", "50", "80")) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  scale_x_continuous(limits = c(-20, 660), breaks = seq(0, 660, 220)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = expression(bold("Biomass: pot volume (g L"^"-1"*")")),
       fill = "Shade cover (%)") +
  guides(linetype = "none", color = "none") +
  pubtheme +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
bvr.cotton


jpeg("../../compile_dissertation/ch2_LxN_Greenhouse/figs/figs1_bvr.jpg",
     height = 7, width = 16, units = 'in', res = 600)
ggarrange(bvr.cotton, bvr.soy,
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "right", 
          align = "hv")
dev.off()

###########################################################
# Root mass fraction plots
########################################################### 
rmf.soy <- ggplot(data = subset(df, spp == "Soybean"), 
                  aes(x = n.ppm,
                      y = rmf,
                      color = factor(shade.cover))) +
  geom_jitter(aes(color = factor(shade.cover)), 
              size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 0.218, yend = 1.33E-05*630 + 0.218), 
               color = "#E69F00", size = 3, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 630, y = 0.168, yend = -5.02E-06*630 + 0.168), 
               color = "#56B4E9", size = 3, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 630, y = 0.142, yend = -1.72E-05*630 + 0.142), 
               color = "#009E73", size = 3, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 630, y = 0.110, yend = -3.56E-05*630 + 0.110), 
               color = "#CC79A7", size = 3, linetype = "dashed") +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = NULL,
       color = "Shade cover (%)") +
  pubtheme +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA))) +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
rmf.soy

rmf.cotton <- ggplot(data = subset(df, spp == "Cotton"), 
                     aes(x = n.ppm,
                         y = rmf,
                         color = factor(shade.cover))) +
  geom_jitter(aes(color = factor(shade.cover)), 
              size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 0.190, yend = -5.70E-05*630 + 0.190), 
               color = "#E69F00", size = 3, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 630, y = 0.159, yend = -4.88E-05*630 + 0.159), 
               color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.141, yend = -4.40E-05*630 + 0.141), 
               color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.118, yend = -3.77E-05*630 + 0.118), 
               color = "#CC79A7", size = 3, linetype = "dashed") +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1)) +
  labs(x = expression(bold("Nitrogen fertilization (ppm)")),
       y = expression(bold("Root mass fraction (unitless)")),
       color = "Shade cover (%)") +
  pubtheme +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA))) +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
rmf.cotton

r <- ggarrange(rmf.cotton, rmf.soy, nrow = 1, ncol = 2,
               common.legend = TRUE,
               legend = "right", align = "hv", 
               labels = NULL, 
               font.label = list(size = 18))
r
ggsave("/Users/eaperkowski/git/manuscripts/N_demand_paper/figs/figS2.rmf.png",
       r,
       width = 16,
       height = 7,
       dpi = "retina")

###########################################################
# Whole plant biomass plots
########################################################### 
tot.soy <- ggplot(data = subset(df, spp == "Soybean"), 
                  aes(x = n.ppm,
                      y = total.bio,
                      color = factor(shade.cover))) +
  geom_jitter(aes(color = factor(shade.cover)), 
              size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 2.142, yend = 2.42E-03*630 + 2.142), 
               color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 1.734, yend = 1.71E-03*630 + 1.734), 
               color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 1.486, yend = 1.29E-03*630 + 1.486), 
               color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 1.150, yend = 7.55E-04*630 + 1.150), 
               color = "#CC79A7", size = 3) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 2)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = NULL,
       color = "Shade cover (%)") +
  pubtheme +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA))) +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
tot.soy

tot.cotton <- ggplot(data = subset(df, spp == "Cotton"), 
                     aes(x = n.ppm,
                         y = total.bio,
                         color = factor(shade.cover))) +
  geom_jitter(aes(color = factor(shade.cover)), 
              size = 4.5, alpha = 0.5, width = 23.625) +
  geom_segment(aes(x = 0, xend = 630, y = 2.373, yend = 2.38E-03*630 + 2.373), 
               color = "#E69F00", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 1.698, yend = 1.11E-03*630 + 1.698), 
               color = "#56B4E9", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 1.358, yend = 6.00E-04*630 + 1.358), 
               color = "#009E73", size = 3) +
  geom_segment(aes(x = 0, xend = 630, y = 0.971, yend = 1.49E-04*630 + 0.971), 
               color = "#CC79A7", size = 3, linetype = "dashed") +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30",
                                "50", "80")) +
  scale_x_continuous(limits = c(-30, 660), breaks = seq(0, 660, 220)) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 2)) +
  labs(x = expression(bold("Nitrogen fertilization (ppm)")),
       y = expression(bold("Total plant biomass (g)")),
       color = "Shade cover (%)") +
  pubtheme +
  guides(fill = FALSE,
         linetype = FALSE,
         color = guide_legend(override.aes = list(fill = NA))) +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
tot.cotton

s <- ggarrange(tot.cotton, tot.soy, nrow = 1, ncol = 2,
               common.legend = TRUE,
               legend = "right", align = "hv", 
               labels = NULL, 
               font.label = list(size = 18))
s
ggsave("/Users/eaperkowski/git/manuscripts/N_demand_paper/figs/figS3.totalBiomass.png",
       s,
       width = 16,
       height = 7,
       dpi = "retina")
