###########################################################
## Load libraries
###########################################################
library(tidyverse)

###########################################################
## Load libraries, create species labels, and convert fixed
## effects to factor
###########################################################
source("https://raw.githubusercontent.com/eaperkowski/LxN_Greenhouse/main/R_files/create_ndemand_metrics.R")

species.label <- c("G. max", "G. hirsutum")
names(species.label) <- c("Soybean", "Cotton")

df$spp <- as.factor(df$spp)
df$n.ppm <- as.factor(df$n.ppm)
df$shade.cover <- as.factor(df$shade.cover)

###########################################################
## Add species mean values
###########################################################
# Add species mean columns
df <- df %>% group_by(spp) %>%
  mutate(ncost.mean = mean(n.cost, na.rm = TRUE),
         n.acq.mean = mean(n.acq, na.rm = TRUE),
         root.carbon.mean = mean(root.carbon.mass, na.rm = TRUE))

###########################################################
## Add theme
###########################################################
theme <- theme(panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(size = 14,
                                         face = "italic"),
               panel.border = element_rect(size = 2, fill = NA),
               axis.text = element_text(size = 14),
               axis.title = element_text(size = 14,
                                         face = "bold"),
               legend.box.background = element_rect(fill = NA, size = 1),
               legend.key = element_rect(fill = NA),
               legend.text = element_text(size = 12),
               legend.title = element_text(size = 13),
               axis.ticks.length = unit(0.25, "cm"))

###########################################################
## Figure 1: Carbon costs to acquire nitrogen (gC gN-1)
###########################################################
ggplot(df, aes(x = shade.cover,
                 y = n.cost,
                 fill = n.ppm)) +
  geom_hline(aes(yintercept = ncost.mean, group = spp), 
             linetype = "solid", alpha = 0.8,
             color = "red", size = 1) +
  geom_boxplot() +
  scale_fill_grey(start = 0.4, end = 1.0,
                  labels = c("0 ppm", "70 ppm",
                             "210 ppm", "630 ppm")) +
  labs(x = "Shade cover (%)",
       y = expression(bold("Carbon cost to acquire nitrogen (gC gN"^-1~")")),
       fill = "Nitrogen fertilization") +
  theme +
  facet_grid(~ spp, labeller = labeller(spp = species.label))

###########################################################
## Figure 2: Whole-plant nitrogen mass (gN)
###########################################################
ggplot(df, aes(x = shade.cover,
                 y = n.acq,
                 fill = n.ppm)) +
  geom_hline(aes(yintercept = n.acq.mean, group = spp), 
             linetype = "solid", alpha = 0.9,
             color = "red", size = 1) +
  geom_boxplot() +
  scale_fill_grey(start = 0.4, end = 1.0,
                  labels = c("0 ppm", "70 ppm",
                             "210 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 0.16),
                     breaks = seq(0, 0.16, 0.04)) +
  scale_x_discrete(limits = c("80", "50", "30", "0")) +
  labs(x = "Shade cover (%)",
       y = "Whole plant nitrogen mass (gN)",
       fill = "Nitrogen fertilization") +
  theme +
  facet_grid(~ spp, labeller = labeller(spp = species.label))

###########################################################
## Figure 3: Root carbon mass (gC)
###########################################################
ggplot(df, aes(x = shade.cover,
                 y = root.carbon.mass,
                 fill = n.ppm)) +
  geom_hline(aes(yintercept = root.carbon.mean, group = spp),
             linetype = "solid", alpha = 0.9,
             color = "red", size = 1) +
  geom_boxplot() +
  scale_fill_grey(start = 0.4, end = 1.0,
                  labels = c("0 ppm", "70 ppm",
                             "210 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 0.4),
                     breaks = seq(0, 0.4, 0.1)) +
  scale_x_discrete(limits = c("80", "50", "30", "0")) +
  labs(x = "Shade cover (%)",
       y = "Root carbon mass (gC)",
       fill = "Nitrogen fertilization") +
  theme +
  facet_grid(~ spp, labeller = labeller(spp = species.label))

###########################################################
## Figure 4: Root nodule mass (g)
###########################################################
ggplot(subset(df, spp == "Soybean"),
       aes(x = shade.cover,
           y = nod.wt,
           fill = n.ppm)) +
  geom_boxplot() +  
  scale_fill_grey(start = 0.4, end = 1.0,
                  labels = c("0 ppm", "70 ppm",
                             "210 ppm", "630 ppm")) +
  scale_x_discrete(limits = c("80", "50", "30", "0")) +
  labs(x = "Shade cover (%)",
       y = "Root nodule mass (g)",
       fill = "Nitrogen fertilization") +
  theme

