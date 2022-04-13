# ------------------------------------------------------------------------------
# Read in our individual data
# ------------------------------------------------------------------------------

library(lme4)
library(tidyverse)
all <- readRDS(cp_path("analysis/data/derived/all.rds"))
feed <- readRDS(cp_path("analysis/data/derived/feed.rds"))
all$cut_pfemia_mbf <- cut(all$pfemia_mbf, c(0,1,1000, 10000, 100000), include.lowest = TRUE)
feed$cut_pfemia_mbf <- cut(feed$pfemia_mbf, c(0,1,1000, 10000, 100000), include.lowest = TRUE)
all$thromb <- factor(c("No Thromb.", "Thromb.")[all$thromb+1], levels = c("No Thromb.", "Thromb."))
all$anemia <- factor(c("No Anemia", "Anemia")[all$anemia+1], levels = c("No Anemia", "Anemia"))

# ------------------------------------------------------------------------------
# Adjusted Models for Infection in our individual data
# ------------------------------------------------------------------------------

asex_mod <- lme4::glmer(
  pf ~ 1 + gender + agecateg2 + thromb + anemia + sickler + (1|stringid) + (1|homestead) + (1|cluster),
  data = all,
  family = "binomial", control=glmerControl(optCtrl=list(maxfun=2e4))
)
asex_mod <- update(asex_mod,start=getME(sex_mod,c("theta","fixef")),control=glmerControl(optCtrl=list(maxfun=1e5)))

sex_mod <- lme4::glmer(
  gametocyte ~ 1 + gender + agecateg2 + thromb + anemia + sickler + cut_pfemia_mbf + (1|stringid) + (1|homestead) + (1|cluster),
  data = all,
  family = "binomial", control=glmerControl(optCtrl=list(maxfun=2e4))
)
sex_mod <- update(sex_mod,start=getME(sex_mod,c("theta","fixef")),control=glmerControl(optCtrl=list(maxfun=2e4)))

# ------------------------------------------------------------------------------
# aOR Plots
# ------------------------------------------------------------------------------

nms <- setNames(
  c("Male","6-11 years","12-24 years",">24 years","Thrombocytopenia",
    "Anemia","Sickle = AS","Sickle = SS","Sickle = Undetermined"),
  rownames(summary(asex_mod)$coefficients)[-1]
)
asex_or <- or_plot(asex_mod, 50, c(30,70), c(0.5,1,1.5,2), nms, order = c(1:9)) +
  scale_x_log10(limits = c(0.01, 120), breaks = c(0.1, 1, 10)) +
  xlab("Adjusted Odds Ratio (aOR) for Asexual Stage Infection")
save_figs("asexual_ors", cowplot::plot_grid(asex_or, labels = "a"), width = 8, height = 6)

nms <- setNames(
  c("Male","6-11 years","12-24 years",">24 years","Thrombocytopenia",
    "Anemia","Sickle = AS","Sickle = SS","Sickle = Undetermined",
    "Parasitemia: 1-1,000p/uL","Parasitemia: 1,000-10,000p/uL","Parasitemia: >10,000p/uL"),
  rownames(summary(sex_mod)$coefficients)[-1]
)
sex_or <- or_plot(sex_mod, 2000, c(1400,5000), c(0.5,1,1.5,2), nms, order = c(1, 2:4, 10:12, 5, 6, 7:9)) +
  scale_x_log10(limits = c(0.1, 10000), breaks = c(0.1, 1, 10, 100))  +
  xlab("Adjusted Odds Ratio (aOR) for Sexual Stage Infection")
save_figs("sexual_ors", cowplot::plot_grid(sex_or, labels = "b"), width = 8, height = 6)


