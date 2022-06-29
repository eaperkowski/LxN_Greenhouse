###########################################################
# Load libraries
###########################################################
library(lme4)
library(car)
library(emmeans)
library(tidyverse)

emm_options(opt.digits = FALSE)

###########################################################
# Load data frame produced by "create_ndemand_metrics.R"
###########################################################
source("LxN_ncost_create_metrics.R")

###########################################################
# Remove rows with missing biomass or nitrogen data
###########################################################
df <- df %>%
        filter(complete.cases(stem.wt, leaves.wt, n.stem))

###########################################################
# Mixed effect model structure
###########################################################
#   Fixed  effects:
#     - shade.cover   - continuous
#     - n.ppm         - continuous
#
#   Random effects:
#     - block         - three blocks

#   Response variables:
#     - n.acq         - whole plant nitrogen mass
#     - root.carbon   - total root carbon mass
#     - n.cost        - root carbon mass / whole plant nitrogen mass
#     - nod.wgt       - root nodule weight

###########################################################
# Carbon cost to acquire nitrogen (g C g-1 N) - G. max
###########################################################
df$n.cost[c(316, 332)] <- NA

ncost.soy <- lmer(sqrt(n.cost) ~ shade.cover * n.ppm + (1 | block),
                  data = subset(df, spp == "Soybean"))

# Check normality assumptions
plot(ncost.soy)
qqnorm(residuals(ncost.soy))
qqline(residuals(ncost.soy))
shapiro.test(residuals(ncost.soy))
outlierTest(ncost.soy)

# Model results
summary(ncost.soy)
Anova(ncost.soy)

# Post-hoc analyses
## For slope
test(emtrends(ncost.soy, 
              ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list()))

## For intercept
emmeans(ncost.soy, 
        ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)))

## For back-transformed slope
test(emtrends(ncost.soy, 
              ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(),
              transform = "response"))

## For back-transformed intercept
emmeans(ncost.soy, 
        ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        type = "response")

###########################################################
# Carbon cost to acquire nitrogen (g C g-1 N) - G. hirsutum
###########################################################
df$n.cost[c(44, 271, 275, 331)] <- NA

ncost.cotton <- lmer(log(n.cost) ~ shade.cover * n.ppm + (1 | block),
                     data = subset(df, spp == "Cotton"))

# Check normality assumptions
plot(ncost.cotton)
qqnorm(residuals(ncost.cotton))
qqline(residuals(ncost.cotton))
shapiro.test(residuals(ncost.cotton))
outlierTest(ncost.cotton)

# Model results
summary(ncost.cotton)
Anova(ncost.cotton)

# Pairwise comparisons
## For slope
test(emtrends(ncost.cotton, 
              ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list()))

## For intercept
emmeans(ncost.cotton, 
        ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)))

## For back-transformed slope
test(emtrends(ncost.cotton, 
              ~shade.cover, 
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(),
              transform = "response"))

## For back-transformed intercept
emmeans(ncost.cotton, 
        ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        type = "response")

###########################################################
# Whole plant nitrogen mass (g N) - G. max
###########################################################
n.acq.soy <- lmer(sqrt(n.acq) ~ shade.cover * n.ppm + (1 | block), 
                  data = subset(df, spp == "Soybean"))

# Check normality assumptions
plot(n.acq.soy)
qqnorm(residuals(n.acq.soy))
qqline(residuals(n.acq.soy))
shapiro.test(residuals(n.acq.soy))
outlierTest(n.acq.soy)

# Model results
summary(n.acq.soy)
Anova(n.acq.soy)

# Pairwise comparisons
## For slope
test(emtrends(n.acq.soy, 
              ~shade.cover,
              var = "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80)), 
              options = list()))

## For intercept
emmeans(n.acq.soy, 
        ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)))

## For back=transformed slope
test(emtrends(n.acq.soy, 
              ~shade.cover,
              var = "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80)), 
              options = list()))

## For back-transformed intercept
emmeans(n.acq.soy, 
        ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        type = "response")

###########################################################
# Whole plant nitrogen mass (g N) - G. hirsutum
###########################################################
df$n.acq[c(286)] <- NA

n.acq.cotton <- lmer(log(n.acq) ~ shade.cover * n.ppm + (1 | block), 
                     data = subset(df, spp == "Cotton"))

# Check normality assumptions
plot(n.acq.cotton)
qqnorm(residuals(n.acq.cotton))
qqline(residuals(n.acq.cotton))
shapiro.test(residuals(n.acq.cotton))
outlierTest(n.acq.cotton)

# Model results
summary(n.acq.cotton)
Anova(n.acq.cotton)

# Pairwise comparisons
## For slope
test(emtrends(n.acq.cotton, 
              ~shade.cover,
              var = "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80)), 
              options = list()))

## For intercept
emmeans(n.acq.cotton, 
        ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)))

## For back-transformed slope
test(emtrends(n.acq.cotton, 
              ~ shade.cover,
              var = "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80)), 
              options = list(),
              transform = "response"))

## For back-transformed intercept
emmeans(n.acq.cotton, 
        ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        type = "response")

###########################################################
# Total root carbon (g C) - G. max
###########################################################
df$root.carbon.mass[c(86)] <- NA

root.carbon.soy <- lmer(sqrt(root.carbon.mass) ~ shade.cover * n.ppm +
                                (1 | block), 
                        data = subset(df, spp == "Soybean"))

# Check normality assumptions
plot(root.carbon.soy)
qqnorm(residuals(root.carbon.soy))
qqline(residuals(root.carbon.soy))
shapiro.test(residuals(root.carbon.soy))
outlierTest(root.carbon.soy)

# Model results
summary(root.carbon.soy)
Anova(root.carbon.soy)

# Pairwise comparisons
## For slope
test(emtrends(root.carbon.soy, 
              ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list()))

## For intercept
emmeans(root.carbon.soy, 
        ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        options = list())

## For back-transformed slope
test(emtrends(root.carbon.soy, 
              ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(),
              transform = "response"))

## For back-transformed intercept
emmeans(root.carbon.soy, 
        ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        options = list(),
        type = "response")

###########################################################
# Total root carbon (g C) - G. hirsutum
###########################################################
df$root.carbon.mass[44] <- NA

root.carbon.cotton <- lmer(sqrt(root.carbon.mass) ~ shade.cover * n.ppm +
                                   (1 | block), 
                           data = subset(df, spp == "Cotton"))

# Check normality assumptions
plot(root.carbon.cotton)
qqnorm(residuals(root.carbon.cotton))
qqline(residuals(root.carbon.cotton))
shapiro.test(residuals(root.carbon.cotton))
outlierTest(root.carbon.cotton)

# Model results
summary(root.carbon.cotton)
Anova(root.carbon.cotton)

# Pairwise comparisons
## For slope
test(emtrends(root.carbon.cotton, 
              ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list()))

## For intercept
emmeans(root.carbon.cotton, 
        ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        options = list())

## For back-transformed slope
test(emtrends(root.carbon.cotton, 
              ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(),
              transform = "response"))

## For back-transformed intercept
emmeans(root.carbon.cotton, ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        options = list(),
        transform = "response")

###########################################################
# Root nodule weight (g)
###########################################################
df$nod.wt[c(10, 203)] <- NA

nod.wgt <- lmer(sqrt(nod.wt) ~ shade.cover * n.ppm + (1 | block), 
                data = subset(df, spp == "Soybean" & nod.wt > 0 ))

# Check normality assumptions
plot(nod.wgt)
qqnorm(residuals(nod.wgt))
qqline(residuals(nod.wgt))
shapiro.test(residuals(nod.wgt))
outlierTest(nod.wgt)

# Model results
summary(nod.wgt)
Anova(nod.wgt)

# Pairwise comparisons
## For slope
test(emtrends(nod.wgt, 
              ~shade.cover, 
              var = "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80)), 
              options = list()))

## For intercept
emmeans(nod.wgt, 
        ~shade.cover,
        var = "n.ppm", 
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        options = list())

## For back-transformed slope
test(emtrends(nod.wgt, 
              ~shade.cover, 
              var = "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80)), 
              options = list(),
              transform = "response"))

## For back-transformed intercept
emmeans(nod.wgt, 
        ~shade.cover,
        var = "n.ppm", 
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        options = list(),
        type = "response")

###########################################################
# Root nodule : root biomass ratio
###########################################################
nod.root <- lmer(sqrt(nod.root.ratio) ~ shade.cover * n.ppm + (1 | block),
                 data = subset(df, spp == "Soybean"))

# Check normality assumptions
plot(nod.root)
qqnorm(residuals(nod.root))
qqline(residuals(nod.root))
shapiro.test(residuals(nod.root))
outlierTest(nod.root)

# Model results
summary(nod.root)
Anova(nod.root)

# Pairwise comparisons
## For slope
test(emtrends(nod.root, 
              ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list()))

## For intercept
emmeans(nod.root, 
        ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                    shade.cover = c(0, 30, 50, 80)),
        options = list())

## For back-transformed slope
test(emtrends(nod.root, 
              ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(),
              transform = "response"))

## For back-transformed intercept
emmeans(nod.root, 
        ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        options = list(),
        type = "response")
