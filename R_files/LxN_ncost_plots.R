###########################################################
## Load libraries
###########################################################
library(tidyverse)
library(ggpubr)
library(emmeans)
library(lme4)

###########################################################
## Load dataset, create species labels, and colorblind
## friendly palette
###########################################################
source("LxN_ncost_create_metrics.R")
## Source path assumes `LxN_ncost_plots.R` is in root directory

species.label <- c("G. max", "G. hirsutum")
names(species.label) <- c("Soybean", "Cotton")
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

###########################################################
## Remove outliers per analyses (Bonferroni p<0.05)
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
# Remove rows with missing biomass or nitrogen data
###########################################################
df <- df %>%
  filter(complete.cases(stem.wt, leaves.wt, n.stem))

###########################################################
## Carbon cost to acquire nitrogen
###########################################################
df$n.cost[c(316, 332)] <- NA
ncost.soy <- lmer(sqrt(n.cost) ~ shade.cover * n.ppm + (1 | block),
                  data = subset(df, spp == "Soybean"))
test(emtrends(ncost.soy, ~shade.cover, "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80))))

df$n.cost[c(44, 271, 275, 331)] <- NA
ncost.cotton <- lmer(log(n.cost) ~ shade.cover * n.ppm + (1 | block),
                     data = subset(df, spp == "Cotton"))
test(emtrends(ncost.cotton, ~shade.cover, "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80))))

## Emmean fxns for regression lines + error ribbons
ncost.soy.regline <- data.frame(emmeans(ncost.soy, ~shade.cover, "n.ppm",
                                    at = list(n.ppm = seq(0, 630, 5),
                                              shade.cover = c(0, 30, 50, 80)),
                                    type = "response")) %>%
  mutate(linetype = ifelse(shade.cover == 0, "dashed", "solid"))
ncost.cotton.regline <- data.frame(emmeans(ncost.cotton, ~shade.cover, "n.ppm",
                                        at = list(n.ppm = seq(0, 630, 5),
                                                  shade.cover = c(0, 30, 50, 80)),
                                        type = "response"))

ncost.nppm.cotton <- ggplot(data = subset(df, spp == "Cotton"), 
                     aes(x = n.ppm, y = n.cost)) +
  geom_jitter(aes(fill = factor(shade.cover)), 
              size = 4, alpha = 0.75, width = 20, shape = 21) +
  geom_smooth(data = ncost.cotton.regline,
              aes(y = response, color = factor(shade.cover)), 
              size = 4, se = FALSE) +
  #geom_ribbon(data = ncost.cotton.regline,
  #            aes(y = response, ymin = lower.CL, ymax = upper.CL, 
  #                fill = factor(shade.cover)), 
  #            size = 1.5, alpha = 0.25) +
  scale_color_manual(values = cbbPalette,
                    labels = c("0", "30", "50", "80")) +
  scale_fill_manual(values = cbbPalette,
                     labels = c("0", "30", "50", "80")) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  scale_x_continuous(limits = c(-20, 660), breaks = seq(0, 660, 220)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = expression(bold("Carbon cost to acquire nitrogen (gC gN"^"-1"*")")),
       fill = "Shade cover (%)") +
  guides(color = "none") +
  pubtheme +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
ncost.nppm.cotton

ncost.nppm.soy <- ggplot(data = subset(df, spp == "Soybean"), 
                         aes(x = n.ppm,y = n.cost)) +
  geom_jitter(aes(fill = factor(shade.cover)), 
              size = 4, alpha = 0.75, width = 20, shape = 21) +
  geom_smooth(data = ncost.soy.regline,
              aes(y = response, color = factor(shade.cover),
                  linetype = linetype), 
              size = 4, se = FALSE) +
  #geom_ribbon(data = ncost.soy.regline,
  #            aes(y = response, ymin = lower.CL, ymax = upper.CL, 
  #                fill = factor(shade.cover)), 
  #            size = 1.5, alpha = 0.25) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30", "50", "80")) +
  scale_fill_manual(values = cbbPalette,
                     labels = c("0", "30", "50", "80")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  scale_x_continuous(limits = c(-20, 660), breaks = seq(0, 660, 220)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = NULL,
       fill = "Shade cover (%)") +
  guides(linetype = "none", color = "none") +
  pubtheme +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
ncost.nppm.soy

png("../corrigendum_figs/fig1_ncost.png",
    height = 7, width = 16, units = "in", res = 600)
ggarrange(ncost.nppm.cotton, ncost.nppm.soy, common.legend = TRUE, legend = "right",
          align = "hv", font.label = list(size = 18))
dev.off()

###########################################################
## Whole-plant nitrogen mass
###########################################################
n.acq.soy <- lmer(sqrt(n.acq) ~ shade.cover * n.ppm + (1 | block), 
                  data = subset(df, spp == "Soybean"))
test(emtrends(n.acq.soy, ~shade.cover, "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80))))

df$n.acq[c(286)] <- NA
n.acq.cotton <- lmer(log(n.acq) ~ shade.cover * n.ppm + (1 | block), 
                     data = subset(df, spp == "Cotton"))
test(emtrends(n.acq.cotton, ~shade.cover, "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80))))


## Emmean fxns for regression lines + error ribbons
nacq.soy.regline <- data.frame(emmeans(n.acq.soy, ~shade.cover, "n.ppm",
                                        at = list(n.ppm = seq(0, 630, 5),
                                                  shade.cover = c(0, 30, 50, 80)),
                                        type = "response"))
nacq.cotton.regline <- data.frame(emmeans(n.acq.cotton, ~shade.cover, "n.ppm",
                                           at = list(n.ppm = seq(0, 630, 5),
                                                     shade.cover = c(0, 30, 50, 80)),
                                           type = "response"))


nacq.nppm.cotton <- ggplot(data = subset(df, spp == "Cotton"), 
                           aes(x = n.ppm, y = n.acq)) +
  geom_jitter(aes(fill = factor(shade.cover)), 
              size = 4, alpha = 0.75, width = 20, shape = 21) +
  geom_smooth(data = nacq.cotton.regline,
              aes(y = response, color = factor(shade.cover)), 
              size = 4, se = FALSE) +
  #geom_ribbon(data = ncost.soy.regline,
  #            aes(y = response, ymin = lower.CL, ymax = upper.CL, 
  #                fill = factor(shade.cover)), 
  #            size = 1.5, alpha = 0.25) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30", "50", "80")) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("0", "30", "50", "80")) +
  scale_y_continuous(limits = c(0, 0.16), breaks = seq(0, 0.16, 0.04)) +
  scale_x_continuous(limits = c(-20, 660), breaks = seq(0, 660, 220)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = "Whole plant nitrogen biomass (g N)",
       fill = "Shade cover (%)") +
  pubtheme +
  guides(color = "none") +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
nacq.nppm.cotton

nacq.nppm.soy <- ggplot(data = subset(df, spp == "Soybean"), 
                       aes(x = n.ppm, y = n.acq)) +
  geom_jitter(aes(fill = factor(shade.cover)), 
              size = 4, alpha = 0.75, width = 20, shape = 21) +
  geom_smooth(data = nacq.soy.regline,
              aes(y = response, color = factor(shade.cover)), 
              size = 4, se = FALSE) +
  #geom_ribbon(data = nacq.soy.regline,
  #            aes(y = response, ymin = lower.CL, ymax = upper.CL, 
  #                fill = factor(shade.cover)), 
  #            size = 1.5, alpha = 0.25) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30", "50", "80")) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("0", "30", "50", "80")) +
  scale_y_continuous(limits = c(0, 0.16), breaks = seq(0, 0.16, 0.04)) +
  scale_x_continuous(limits = c(-20, 660), breaks = seq(0, 660, 220)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = NULL,
       fill = "Shade cover (%)") +
  guides(color = "none") +
  pubtheme +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
nacq.nppm.soy

png("../corrigendum_figs/fig2_nacq.png",
    height = 7, width = 16, units = "in", res = 600)
ggarrange(nacq.nppm.cotton, nacq.nppm.soy, common.legend = TRUE, 
          legend = "right", align = "hv", font.label = list(size = 18))
dev.off()

###########################################################
## Root carbon mass
###########################################################
df$root.carbon.mass[c(86)] <- NA
cbg.soy <- lmer(sqrt(root.carbon.mass) ~ shade.cover * n.ppm + (1 | block), 
                        data = subset(df, spp == "Soybean"))
test(emtrends(cbg.soy, ~shade.cover, "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80))))

df$root.carbon.mass[44] <- NA
cbg.cotton <- lmer(sqrt(root.carbon.mass) ~ shade.cover * n.ppm + (1 | block), 
                           data = subset(df, spp == "Cotton"))
test(emtrends(cbg.cotton, ~shade.cover, "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80))))

## Emmean fxns for regression lines + error ribbons
cbg.soy.regline <- data.frame(emmeans(cbg.soy, ~shade.cover, "n.ppm",
                                        at = list(n.ppm = seq(0, 630, 5),
                                                  shade.cover = c(0, 30, 50, 80)),
                                        type = "response")) %>%
  mutate(linetype = ifelse(shade.cover == 80, "dashed", "solid"))

cbg.cotton.regline <- data.frame(emmeans(cbg.cotton, ~shade.cover, "n.ppm",
                                           at = list(n.ppm = seq(0, 630, 5),
                                                     shade.cover = c(0, 30, 50, 80)),
                                           type = "response")) %>%
  mutate(linetype = ifelse(shade.cover == 50 | shade.cover == 80, "dashed", "solid"))


cbg.nppm.cotton <- ggplot(data = subset(df, spp == "Cotton"), 
                          aes(x = n.ppm, y = root.carbon.mass)) +
  geom_jitter(aes(fill = factor(shade.cover)), 
              size = 4, alpha = 0.75, width = 20, shape = 21) +
  geom_smooth(data = cbg.cotton.regline,
              aes(y = response, color = factor(shade.cover), linetype = linetype), 
              size = 4, se = FALSE) +
  #geom_ribbon(data = ncost.soy.regline,
  #            aes(y = response, ymin = lower.CL, ymax = upper.CL, 
  #                fill = factor(shade.cover)), 
  #            size = 1.5, alpha = 0.25) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30", "50", "80")) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("0", "30", "50", "80")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2)) +
  scale_x_continuous(limits = c(-20, 660), breaks = seq(0, 660, 220)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = expression(bold("Root carbon biomass (g C)")),
       color = "Shade cover (%)", fill = "Shade cover (%)") +
  guides(linetype = "none", color = "none") +
  pubtheme +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
cbg.nppm.cotton

cbg.nppm.soy <- ggplot(data = subset(df, spp == "Soybean"), 
                              aes(x = n.ppm, y = root.carbon.mass)) +
  geom_jitter(aes(fill = factor(shade.cover)), 
              size = 4, alpha = 0.75, width = 20, shape = 21) +
  geom_smooth(data = cbg.soy.regline,
              aes(y = response, color = factor(shade.cover), linetype = linetype), 
              size = 4, se = FALSE) +
  #geom_ribbon(data = ncost.soy.regline,
  #            aes(y = response, ymin = lower.CL, ymax = upper.CL, 
  #                fill = factor(shade.cover)), 
  #            size = 1.5, alpha = 0.25) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30", "50", "80")) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("0", "30", "50", "80")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2)) +
  scale_x_continuous(limits = c(-20, 660), breaks = seq(0, 660, 220)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = NULL,
       color = "Shade cover (%)", fill = "Shade cover (%)") +
  guides(linetype = "none", color = "none") +
  pubtheme  +
  facet_grid(.~spp, labeller = labeller(spp = species.label))
  
cbg.nppm.soy

png("../corrigendum_figs/fig3_rootCarbon.png",
    height = 7, width = 16, units = "in", res = 600)
ggarrange(cbg.nppm.cotton, cbg.nppm.soy, common.legend = TRUE, 
          legend = "right", align = "hv", font.label = list(size = 18))
dev.off()

###########################################################
## Root nodule weight
###########################################################
df$nod.wt[c(10, 203)] <- NA

nod.wgt <- lmer(sqrt(nod.wt) ~ shade.cover * n.ppm + (1 | block), 
                data = subset(df, spp == "Soybean" & nod.wt > 0))
test(emtrends(nod.wgt, ~shade.cover, "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80))))

## Emmean fxns for regression lines + error ribbons
nod.wgt.regline <- data.frame(emmeans(nod.wgt, ~shade.cover, "n.ppm",
                                        at = list(n.ppm = seq(0, 630, 5),
                                                  shade.cover = c(0, 30, 50, 80)),
                                        type = "response"))

nod.weight.ppm <- ggplot(data = subset(df, spp == "Soybean" & nod.wt > 0), 
                           aes(x = n.ppm, y = nod.wt)) +
  geom_jitter(aes(fill = factor(shade.cover)), 
              size = 4, alpha = 0.75, width = 20, shape = 21) +
  geom_smooth(data = nod.wgt.regline,
              aes(y = response, color = factor(shade.cover)), 
              size = 4, se = FALSE) +
  scale_color_manual(values = cbbPalette,
                     labels = c("0", "30", "50", "80")) +
  scale_fill_manual(values = cbbPalette,
                     labels = c("0", "30", "50", "80")) +
  scale_y_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, 0.05)) +
  scale_x_continuous(limits = c(-20, 660), breaks = seq(0, 660, 220)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = "Root nodule biomass (g)",
       fill = "Shade cover (%)") +
  guides(color = "none") +
  pubtheme
nod.weight.ppm

###########################################################
## Root nodule biomass : root biomass
###########################################################
nod.root <- lmer(sqrt(nod.root.ratio) ~ shade.cover * n.ppm + (1 | block),
                 data = subset(df, spp == "Soybean"))
test(emtrends(nod.root, ~shade.cover, "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80))))

## Emmean fxns for regression lines + error ribbons
nod.root.regline <- data.frame(emmeans(nod.root, ~shade.cover, "n.ppm",
                                      at = list(n.ppm = seq(0, 630, 5),
                                                shade.cover = c(0, 30, 50, 80)),
                                      type = "response"))

nod.root.ppm <- ggplot(data = subset(df, spp == "Soybean" & nod.root.ratio > 0),
                       aes(x = n.ppm, y = nod.root.ratio)) +
  geom_jitter(aes(fill = factor(shade.cover)), 
              size = 4, width = 20, alpha = 0.75, shape = 21) +
  geom_smooth(data = nod.root.regline,
              aes(y = response, color = factor(shade.cover)), 
              size = 4, se = FALSE) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("0", "30", "50", "80")) +
  scale_color_manual(values = cbbPalette,
                    labels = c("0", "30", "50", "80")) +
  scale_x_continuous(limits = c(-20, 660), breaks = seq(0, 660, 220)) +
  labs(x = "Nitrogen fertilization (ppm)",
       y = expression(bold("Nodule: root biomass (g g"^"-1"*")")),
       fill = "Shade cover (%)", color = "Shade cover (%)") +
  pubtheme
nod.root.ppm



png("../corrigendum_figs/fig4_nodwgt.png",
    height = 6.5, width = 16, units = "in", res = 300)
ggarrange(nod.weight.ppm, nod.root.ppm, common.legend = TRUE, legend = "right",
          align = "hv", labels = "AUTO", font.label = list(size = 18))
dev.off()

# png("../../compile_dissertation/ch2_LxN_Greenhouse/figs/fig4_nodwgt.png",
#     height = 4.5, width = 12, units = "in", res = 600)
# ggarrange(nod.weight.ppm, nod.root.ppm, common.legend = TRUE, legend = "right",
#           align = "hv", labels = "AUTO", font.label = list(size = 18))
# dev.off()
