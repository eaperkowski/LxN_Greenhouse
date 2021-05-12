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
  filter(complete.cases(stem.wt, leaves.wt, n.stem))

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
## For bvr-n.ppm slope
emtrends(bvr.cotton, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(), 
              transform = "response")
## For bvr-n.ppm significance level
test(emtrends(bvr.cotton, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(), 
              type = "response"))

## For bvr-n.ppm intercept
emmeans(bvr.cotton, ~n.ppm * shade.cover,
                   at = list(n.ppm = 0,
                             shade.cover = c(0, 30, 50, 80)),
        type = "response")  

###########################################################
# Root mass fraction - G. max
###########################################################
df$rmf[c(86, 90, 316, 332)] <- NA

rmf.soy <- lmer(rmf ~ shade.cover * n.ppm + (1 | block), 
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
## For rmf-n.ppm slope
emtrends(rmf.soy, ~shade.cover,
         var = "n.ppm",
         at = list(shade.cover = c(0, 30, 50, 80)),
         options = list(), 
         transform = "response")
## For rmf-n.ppm significance level
test(emtrends(rmf.soy, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(), 
              type = "response"))

## For rmf-n.ppm intercept
emmeans(rmf.soy, ~n.ppm * shade.cover,
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        type = "response")  

###########################################################
# Root mass fraction - G. hirsutum
###########################################################
df$rmf[c(44, 171, 271, 275, 331)] <- NA

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
## For rmf-n.ppm slope
emtrends(rmf.cotton, ~shade.cover,
         var = "n.ppm",
         at = list(shade.cover = c(0, 30, 50, 80)),
         options = list(), 
         transform = "response")
## For rmf-n.ppm significance level
test(emtrends(rmf.cotton, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(), 
              type = "response"))

## For rmf-n.ppm intercept
emmeans(rmf.cotton, ~n.ppm * shade.cover,
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        type = "response") 

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
## For tot-n.ppm slope
emtrends(tot.soy, ~shade.cover,
         var = "n.ppm",
         at = list(shade.cover = c(0, 30, 50, 80)),
         options = list(), 
         transform = "response")
## For tot-n.ppm significance level
test(emtrends(tot.soy, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(), 
              type = "response"))

## For tot-n.ppm intercept
emmeans(tot.soy, ~n.ppm * shade.cover,
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        type = "response")  

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
## For tot-n.ppm slope
emtrends(tot.cotton, ~shade.cover,
         var = "n.ppm",
         at = list(shade.cover = c(0, 30, 50, 80)),
         options = list(), 
         transform = "response")
## For tot-n.ppm significance level
test(emtrends(tot.cotton, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(), 
              type = "response"))

## For tot-n.ppm intercept
emmeans(tot.cotton, ~n.ppm * shade.cover,
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        type = "response") 