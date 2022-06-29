###########################################################
## Load libraries
###########################################################
library(lme4)
library(car)
library(emmeans)

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
## For bvr:n.ppm slopes & slope significance level
test(emtrends(bvr.soy, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list()))

## For bvr:n.ppm intercepts
emmeans(bvr.soy, ~n.ppm * shade.cover,
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)))

## For back-transformed bvr:n.ppm intercepts
emmeans(bvr.soy, ~n.ppm * shade.cover,
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        type = "response")

## What was the marginal mean for 0% shade:630 ppm N combo?
emmeans(bvr.soy, ~n.ppm*shade.cover,
        at = list(n.ppm = 630,
                  shade.cover = 0))

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
## For bvr-n.ppm slope (figs only)
emtrends(bvr.cotton, ~shade.cover,
         var = "n.ppm",
         at = list(shade.cover = c(0, 30, 50, 80)),
         options = list(),
         transform = "response")

## For bvr-n.ppm slope significance level
test(emtrends(bvr.cotton, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list()))

## For back-transformed bvr-n.ppm intercept
emmeans(bvr.cotton, ~n.ppm * shade.cover,
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        type = "response")

## For bvr-n.ppm intercept
emmeans(bvr.cotton, ~n.ppm * shade.cover,
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)))


## What was the marginal mean for 0% shade:630 ppm N combo?
emmeans(bvr.cotton, ~n.ppm*shade.cover,
        at = list(n.ppm = 630,
                  shade.cover = 0),
        type = "response")
