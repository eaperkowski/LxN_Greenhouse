## Analyses for all tables in supplemental material

###########################################################
## Load libraries
###########################################################
library(lme4)
library(car)
library(emmeans)
library(RVAideMemoire)

###########################################################
## Load data set
###########################################################
source("https://raw.githubusercontent.com/eaperkowski/LxN_Greenhouse/main/R_files/create_ndemand_metrics.R")

###########################################################
## Calculate BVR, RMF, and total biomass
###########################################################
df <- df %>%
  mutate(total.bio = stem.wt + roots.wt + leaves.wt,
         rmf = roots.wt / total.bio,
         bvr = total.bio / 3) %>%
  filter(complete.cases(stem.wt))

###########################################################
# BVR - G. max
###########################################################
bvr.soy <- lmer(bvr ~ shade.cover * n.ppm + (1 | block), 
                data = subset(df, spp == "Soybean"))

# Check normality assumptions
plot(bvr.soy)
qqnorm(residuals(bvr.soy))
qqline(residuals(bvr.soy))
shapiro.test(residuals(bvr.soy))
outlierTest(bvr.soy)

# Model results
summary(bvr.soy)
Anova(bvr.soy)

# Pairwise comparisons
test(emtrends(bvr.soy, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list()))
emmeans(bvr.soy, ~n.ppm * shade.cover,
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)))

###########################################################
# BVR - G. hirsutum
###########################################################
bvr.cotton <- lmer(log(bvr) ~ shade.cover * n.ppm + (1 | block), 
                data = subset(df, spp == "Cotton"))

# Check normality assumptions
plot(bvr.cotton)
qqnorm(residuals(bvr.cotton))
qqline(residuals(bvr.cotton))
shapiro.test(residuals(bvr.cotton))
outlierTest(bvr.cotton)

# Model results
summary(bvr.cotton)
Anova(bvr.cotton)

# Pairwise comparisons
test(emtrends(bvr.cotton, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(), transform = "response"))
bvr.cotton.int <- emmeans(bvr.cotton, ~n.ppm * shade.cover,
                   at = list(n.ppm = 0,
                             shade.cover = c(0, 30, 50, 80)))  
back.emmeans(bvr.cotton.int, transform = "log")

###########################################################
# Root mass fraction - G. max
###########################################################
df$rmf[c(10, 87, 91, 318, 334)] <- NA

rmf.soy <- lmer(log(rmf) ~ shade.cover * n.ppm + (1 | block), 
                data = subset(df, spp == "Soybean"))

# Check normality assumptions
plot(rmf.soy)
qqnorm(residuals(rmf.soy))
qqline(residuals(rmf.soy))
shapiro.test(residuals(rmf.soy))
outlierTest(rmf.soy)

# Model results
summary(rmf.soy)
Anova(rmf.soy)

# Pairwise comparisons
test(emtrends(rmf.soy, ~shade.cover, 
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(),
              transform = "response"))
rmf.soy.int <- emmeans(rmf.soy, ~shade.cover,
                       var = "n.ppm",
                       at = list(n.ppm = 0,
                                 shade.cover = c(0, 30, 50, 80)))
back.emmeans(rmf.soy.int, transform = "log")

###########################################################
# Root mass fraction - G. hirsutum
###########################################################
df$rmf[c(44, 172, 273, 277, 333)] <- NA

rmf.cotton <- lmer(log(rmf) ~ shade.cover * n.ppm + (1 | block), 
                data = subset(df, spp == "Cotton"))

# Check normality assumptions
plot(rmf.cotton)
qqnorm(residuals(rmf.cotton))
qqline(residuals(rmf.cotton))
shapiro.test(residuals(rmf.cotton))
outlierTest(rmf.cotton)

# Model results
summary(rmf.cotton)
Anova(rmf.cotton)

# Pairwise comparisons
test(emtrends(rmf.cotton, ~shade.cover, 
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(),
              transform = "response"))
rmf.cotton.int <- emmeans(rmf.cotton, ~shade.cover,
                          var = "n.ppm",
                          at = list(n.ppm = 0,
                                    shade.cover = c(0, 30, 50, 80)))
back.emmeans(rmf.cotton.int, transform = "log")

###########################################################
# Total biomass - G. max
###########################################################
tot.soy <- lmer(sqrt(total.bio) ~ shade.cover * n.ppm + (1 | block), 
                data = subset(df, spp == "Soybean"))

# Check normality assumptions
plot(tot.soy)
qqnorm(residuals(tot.soy))
qqline(residuals(tot.soy))
shapiro.test(residuals(tot.soy))
outlierTest(tot.soy)

# Model results
summary(tot.soy)
Anova(tot.soy)

# Pairwise comparisons
test(emtrends(tot.soy, ~shade.cover, 
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(),
              transform = "response"))
tot.soy.int <- emmeans(tot.soy, ~shade.cover,
                          var = "n.ppm",
                          at = list(n.ppm = 0,
                                    shade.cover = c(0, 30, 50, 80)))
back.emmeans(tot.soy.int, transform = "sqrt")

###########################################################
# Total  biomass - G. hirsutum
###########################################################
tot.cotton <- lmer(log(total.bio) ~ shade.cover * n.ppm + (1 | block), 
                   data = subset(df, spp == "Cotton"))

# Check normality assumptions
plot(tot.cotton)
qqnorm(residuals(tot.cotton))
qqline(residuals(tot.cotton))
shapiro.test(residuals(tot.cotton))
outlierTest(tot.cotton)

# Model results
summary(tot.cotton)
Anova(tot.cotton)

# Pairwise comparisons
test(emtrends(tot.cotton, ~shade.cover, 
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(),
              transform = "response"))
tot.cotton.int <- emmeans(tot.cotton, ~shade.cover,
                          var = "n.ppm", 
                          at = list(n.ppm = 0,
                                    shade.cover = c(0, 30, 50, 80)))
back.emmeans(tot.cotton.int, transform = "log")