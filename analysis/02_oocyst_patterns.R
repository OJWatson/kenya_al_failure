# ------------------------------------------------------------------------------
# Read in our feed data
# ------------------------------------------------------------------------------

library(tidyverse)
library(lme4)
devtools::load_all()
feed <- readRDS(cp_path("analysis/data/derived/feed.rds"))

############ basic plots of feeding -----------------

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"),...)
}

# check power vs linear (gompertz and hyperbolic also checked and known to be much worse)
# power won't converge right and is with
pow16dfa <- minpack.lm::nlsLM(oocystdens ~ a0 + a1*pf16c^a2, start = list(a0 = 5, a1=1, a2=5), data = feed %>% filter(assay == "DFA"), control = list(maxiter = 100))
lin16dfa <- minpack.lm::nlsLM(oocystdens ~ a0 + a1*pf16c , start = list(a0 = 0.1, a1=0.1), data = feed %>% filter(assay == "DFA"))
which.min(c(AIC(pow16dfa), AIC(lin16dfa)))
dfa16 <- data.frame("assay" = "DFA","gutdissected"=1,"pf16c" = seq(10,40,0.01),
                           "oocystdens" = predict(pow16dfa, type = "response", newdata = data.frame(pf16c = seq(10,40,0.01))))


pow16mfa <- minpack.lm::nlsLM(oocystdens ~ a0 + a1*pf16c^a2, start = list(a0 = 5, a1=1, a2=5), data = feed %>% filter(assay == "MFA"), control = list(maxiter = 100))
lin16mfa <- minpack.lm::nlsLM(oocystdens ~ a0 + a1*pf16c , start = list(a0 = 0.1, a1=0.1), data = feed %>% filter(assay == "MFA"))
which.min(c(AIC(pow16mfa), AIC(lin16mfa)))
mfa16 <- data.frame("assay" = "MFA","gutdissected"=1,"pf16c" = seq(10,40,0.01),
                           "oocystdens" = predict(pow16mfa, type = "response", newdata = data.frame(pf16c = seq(10,40,0.01))))

pow25dfa <- minpack.lm::nlsLM(oocystdens ~ a0 + a1*pf25c^a2, start = list(a0 = 5, a1=1, a2=5), data = feed %>% filter(assay == "DFA"), control = list(maxiter = 100))
lin25dfa <- minpack.lm::nlsLM(oocystdens ~ a0 + a1*pf25c , start = list(a0 = 0.1, a1=0.1), data = feed %>% filter(assay == "DFA"))
which.min(c(AIC(pow25dfa), AIC(lin25dfa)))
dfa25 <- data.frame("assay" = "DFA","gutdissected"=1,"pf25c" = seq(10,40,0.01),
                    "oocystdens" = predict(pow25dfa, type = "response", newdata = data.frame(pf25c = seq(10,40,0.01))))


pow25mfa <- minpack.lm::nlsLM(oocystdens ~ a0 + a1*pf25c^a2, start = list(a0 = 5, a1=1, a2=5), data = feed %>% filter(assay == "MFA"), control = list(maxiter = 100))
lin25mfa <- minpack.lm::nlsLM(oocystdens ~ a0 + a1*pf25c , start = list(a0 = 0.1, a1=0.1), data = feed %>% filter(assay == "MFA"))
which.min(c(AIC(pow25mfa), AIC(lin25mfa)))
mfa25 <- data.frame("assay" = "MFA","gutdissected"=1,"pf25c" = seq(10,40,0.01),
                    "oocystdens" = predict(pow25mfa, type = "response", newdata = data.frame(pf25c = seq(10,40,0.01))))

# pow <- minpack.lm::nlsLM(oocystdens ~ a0 + a1*pf25c^a2 , start = list(a0 = 1, a1=100, a2=-2), data = feed)
# lin <- minpack.lm::nlsLM(oocystdens ~ a0 + a1*pf25c , start = list(a0 = 0.1, a1=0.1), data = feed)
# hyp <- minpack.lm::nlsLM(oocystdens ~ a0 + (a1*pf25c)/(1+a3*pf25c), start = list(a0 = 5, a1=100, a3 = 100), data = feed, lower = c(-Inf, 0, 0))
# gomp <- minpack.lm::nlsLM(oocystdens ~ a0 + (a1*pf25c^a2)/(1+a3*pf25c^a2), start = list(a0 = 10, a1=10, a2=2, a3 = 10), data = feed, lower = c(-Inf, -Inf, 1, -Inf))
#
# which.min(c(AIC(pow), AIC(lin), AIC(hyp), AIC(gomp)))


int_prev <- feed %>%
  filter(oocystprevalence > 0) %>%
  ggplot(aes(oocystdens, oocystprevalence/100, weight = gutdissected, color = assay)) +
  geom_point(alpha = 0.8) +
  binomial_smooth(se = FALSE, lwd = 1.5) +
  coord_trans(x="log10") +
  scale_x_continuous(breaks = c(1,10,100)) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  xlab("Ooyst Density") +
  ylab("Oocyst Prevalence") +
  scale_color_viridis_d(name = "Assay", end = 0.8) +
  theme(panel.grid.minor.x = element_blank())

pf16 <- feed %>%
  filter(!is.na(pf16c)) %>%
  ggplot(aes(pf16c, oocystdens, weight = gutdissected, color = assay)) +
  geom_point(alpha = 0.8) +
  geom_smooth(se = FALSE, method = "lm", lwd = 1.5) +
  theme_bw() +
  xlab("Pfs16 Ct") +
  ylab("Oocyst Density") +
  scale_y_log10() +
  scale_color_viridis_d(name = "Assay", end = 0.8)  +
  theme(panel.grid.minor.x = element_blank())

pf25 <- feed  %>%
  filter(!is.na(pf25c)) %>%
  ggplot(aes(pf25c, oocystdens, weight = gutdissected, color = assay)) +
  geom_point(alpha = 0.8) +
  geom_smooth(se = FALSE, method = "lm", lwd = 1.5) +
  theme_bw() +
  xlab("Pfs25 Ct") +
  ylab("Oocyst Density") +
  scale_color_viridis_d(name = "Assay", end = 0.8) +
  scale_y_log10(limits = c(1, 50)) +
  theme(panel.grid.minor.x = element_blank())


# Multivariate Logistic Regression for odds of mosquito infection
# collapse between 1 and 0 for bet_binomial (https://github.com/glmmTMB/glmmTMB/issues/507#issuecomment-513235407)
feed_mod_dat <- feed %>% mutate(oocystprevalence = replace(oocystprevalence, which(oocystprevalence==100), 99.99),
                                oocystprevalence = replace(oocystprevalence, which(oocystprevalence==0), 0.01))

prev_mod <- glmmTMB::glmmTMB(
  oocystprevalence/100 ~ 1 + visit_prev  + assay + (pf16c) + (pf25c) + subgf + (1|homestead),
  data = feed_mod_dat %>% mutate(assay = factor(assay, levels = c("MFA","DFA")), visit_prev = visit_prev/100),
  family = glmmTMB::beta_family())

nms <- setNames(
  c("Monthly Parasite\nPrevalence","Direct\nFeeding\nAssay", "Pfs16 Ct", "Pfs25 Ct", "Submicroscopic\n Gametocytemia"),
  rownames(summary(prev_mod)$coefficients$cond)[-1]
)
prev_or <- or_plot(prev_mod, 4, c(0.5,0.5), c(0.5,1,1.5,2, 2.5, 3, 3.5), nms, 1:5) +
  xlab("Adjusted Odds Ratio (aOR) for Ooocyst Prevalence")


# Bring plots all together
fig2a <- cowplot::plot_grid(pf16 + theme(legend.position = "none", legend.background = element_rect(color = "black")),
                   pf25 +  theme(legend.position = "none", legend.background = element_rect(color = "black")),
                   int_prev +  theme(legend.position = "none", legend.background = element_rect(color = "black")),
                   cowplot::get_legend(pf16),
                   ncol = 4, rel_widths = c(1,1,1,0.4),
                   labels = c("a", "b", "c")) +
  theme(plot.background = element_rect(fill = "white", color = "white"))


fig2 <- cowplot::plot_grid(fig2a, NA, cowplot::plot_grid(NA, prev_or, rel_widths = c(0.01,6)),
                           ncol = 1,
                           rel_heights = c(1,0.1, 1), labels = c("", "", "d")) +
  theme(plot.background = element_rect(fill = "white"))

save_figs("oocyst_patterns", fig2, width = 12, height = 10)


## oocyst gland prev

oocgld <- feed %>% ggplot(aes(oocystprevalence/100, glandprevalence/100, color = assay)) +
  binomial_smooth(alpha = 0.5) +
  geom_point() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  theme_bw() +
  xlab("Oocyst Prevalence") +
  ylab("Gland Prevalence") +
  scale_color_viridis_d(name = "Assay", end = 0.8) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent)

save_figs("gland_prevalence", oocgld, width = 6, height = 4)


# Numbers

glm(
  oocystprevalence ~ glandprevalence,
  feed %>% mutate(oocystprevalence = oocystprevalence/100, glandprevalence = glandprevalence/100) %>%
    filter(assay == "DFA"),
  family = "binomial") %>%
  broom::tidy(conf.int = TRUE) %>%
  mutate(across(c("estimate", "conf.low", "conf.high"), .fns = exp)
  )

glm(
  oocystprevalence ~ glandprevalence,
  feed %>% mutate(oocystprevalence = oocystprevalence/100, glandprevalence = glandprevalence/100) %>%
    filter(assay == "MFA"),
  family = "binomial") %>%
  broom::tidy(conf.int = TRUE) %>%
  mutate(across(c("estimate", "conf.low", "conf.high"), .fns = exp)
  )


