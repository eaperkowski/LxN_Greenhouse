###########################################################
# Load libraries
###########################################################
library(lme4)
library(car)
library(emmeans)
library(MuMIn)
library(ggplot2)

###########################################################
# Load data frame produced by "create_ndemand_metrics.R"
###########################################################
source("https://raw.githubusercontent.com/eaperkowski/LxN_Greenhouse/main/R_files/create_ndemand_metrics.R")

###########################################################
# Mixed effect model structure
###########################################################
#   Fixed  effects:
#     - shade.cover   - continuous
#     - n.ppm         - continuous
#
#   Random effects:
#     - block         - three blocks
#
#   Response variables:
#     - n.acq         - whole plant nitrogen mass
#     - root.carbon   - total root carbon mass
#     - n.cost        - root carbon mass / whole plant nitrogen mass
#     - nod.wgt       - root nodule weight

###########################################################
# Carbon cost to acquire nitrogen (g C g-1 N) - G. max
###########################################################
df$n.cost[c(87, 321)] <- NA

ncost.soy <- lmer(log(n.cost) ~ shade.cover * n.ppm + (1 | block), 
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
r.squaredGLMM(ncost.soy)

# Post-hoc analyses
test(emtrends(ncost.soy, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(),
              transform = "response"))
emmeans(ncost.soy, ~n.ppm * shade.cover,
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        transform = "response")

###########################################################
# Carbon cost to acquire nitrogen (g C g-1 N) - G. hirsutum
###########################################################
df$n.cost[c(44, 275, 279, 337)] <- NA

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
r.squaredGLMM(ncost.cotton)

# Pairwise comparisons
test(emtrends(ncost.cotton, ~shade.cover, 
         var = "n.ppm",
         at = list(shade.cover = c(0, 30, 50, 80)),
         options = list(),
         transform = "response"))

emmeans(ncost.cotton, ~n.ppm * shade.cover,
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        transform = "response")

###########################################################
# Whole plant nitrogen mass (g N) - G. max
###########################################################
df$n.acq[20] <- NA

n.acq.soy <- lmer(log(n.acq) ~ shade.cover * n.ppm + (1 | block), 
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
r.squaredGLMM(n.acq.soy)

# Pairwise comparisons
test(emtrends(n.acq.soy, ~shade.cover,
              var = "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80)), 
              options = list(),
              transform = "response"))
emmeans(n.acq.soy, ~n.ppm * shade.cover,
        at = list(n.ppm = 0,
                  shade.cover = c(0, 30, 50, 80)),
        transform = "response")

###########################################################
# Whole plant nitrogen mass (g N) - G. hirsutum
###########################################################
df$n.acq[290] <- NA

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
r.squaredGLMM(n.acq.cotton)

# Pairwise comparisons
test(emtrends(n.acq.cotton, ~ shade.cover,
              var = "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80)), 
              options = list(),
              transform = "response"))
emmeans(n.acq.cotton, ~n.ppm*shade.cover,
        at = list(n.ppm = c(0),
                  shade.cover = c(0, 30, 50, 80)),
        transform = "response")

###########################################################
# Total root carbon (g C) - G. max
###########################################################
df$root.carbon.mass[87] <- NA

root.carbon.soy <- lmer(sqrt(root.carbon.mass) ~ shade.cover * n.ppm +
                      (1 | block), data = subset(df, spp == "Soybean"))

# Check normality assumptions
plot(root.carbon.soy)
qqnorm(residuals(root.carbon.soy))
qqline(residuals(root.carbon.soy))
shapiro.test(residuals(root.carbon.soy))
outlierTest(root.carbon.soy)

# Model results
summary(root.carbon.soy)
Anova(root.carbon.soy)
r.squaredGLMM(root.carbon.soy)

# Pairwise comparisons
test(emtrends(root.carbon.soy, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(),
              transform = "response"))

emmeans(root.carbon.soy, ~shade.cover,
         var = "n.ppm",
         at = list(n.ppm = 0,
                   shade.cover = c(0, 30, 50, 80)),
         options = list(),
        transform = "response")

###########################################################
# Total root carbon (g C) - G. hirsutum
###########################################################
df$root.carbon.mass[44] <- NA

root.carbon.cotton <- lmer(sqrt(root.carbon.mass) ~ shade.cover * n.ppm +
                      (1 | block), data = subset(df, spp == "Cotton"))

# Check normality assumptions
plot(root.carbon.cotton)
qqnorm(residuals(root.carbon.cotton))
qqline(residuals(root.carbon.cotton))
shapiro.test(residuals(root.carbon.cotton))
outlierTest(root.carbon.cotton)

# Model results
summary(root.carbon.cotton)
Anova(root.carbon.cotton)
r.squaredGLMM(root.carbon.cotton)

# Pairwise comparisons
test(emtrends(root.carbon.cotton, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(),
              transform = "response"))
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
r.squaredGLMM(nod.wgt)

# Pairwise comparisons
test(emtrends(nod.wgt, ~shade.cover, 
              var = "n.ppm", 
              at = list(shade.cover = c(0, 30, 50, 80)), 
              options = list(),
              transform = "response"))
emmeans(nod.wgt, ~shade.cover, 
          var = "n.ppm", 
          at = list(n.ppm = 0,
                    shade.cover = c(0, 30, 50, 80)), 
          options = list(),
        transform = "response")

###########################################################
# Root nodule : root biomass ratio
###########################################################
df$nod.root.ratio <- df$nod.wt / df$roots.wt

nod.root <- lmer(sqrt(nod.root.ratio) ~ shade.cover * n.ppm + (1 | block), 
                data = subset(df, spp == "Soybean" & nod.wt > 0))

# Check normality assumptions
plot(nod.root)
qqnorm(residuals(nod.root))
qqline(residuals(nod.root))
shapiro.test(residuals(nod.root))
outlierTest(nod.root)

# Model results
summary(nod.root)
Anova(nod.root)
r.squaredGLMM(nod.root)

test(emtrends(nod.root, ~shade.cover,
              var = "n.ppm",
              at = list(shade.cover = c(0, 30, 50, 80)),
              options = list(),
              transform = "response"))
emmeans(nod.root, ~shade.cover,
        var = "n.ppm",
        at = list(n.ppm = 0,
                    shade.cover = c(0, 30, 50, 80)),
        options = list(),
        transform = "response")
