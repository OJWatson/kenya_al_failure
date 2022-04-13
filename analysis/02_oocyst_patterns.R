# ------------------------------------------------------------------------------
# Read in our feed data
# ------------------------------------------------------------------------------

library(tidyverse)
devtools::load_all()
feed <- readRDS(cp_path("analysis/data/derived/feed.rds"))

############ basic plots of feeding -----------------

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"),...)
}

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
  xlab("Pf16 Ct") +
  ylab("Oocyst Density") +
  scale_color_viridis_d(name = "Assay", end = 0.8) +
  theme(panel.grid.minor.x = element_blank())

pf25 <- feed  %>%
  filter(!is.na(pf25c)) %>%
  ggplot(aes(pf25c, oocystdens, weight = gutdissected, color = assay)) +
  geom_point(alpha = 0.8) +
  geom_smooth(se = FALSE, method = "lm", lwd = 1.5) +
  theme_bw() +
  xlab("Pf25 Ct") +
  ylab("Oocyst Prevalence") +
  scale_color_viridis_d(name = "Assay", end = 0.8) +
  theme(panel.grid.minor.x = element_blank())


# Multivariate Logistic Regression for odds of mosquito infection
prev_mod <- lme4::glmer(
  oocystprevalence/100 ~ 1 + visit_prev + (oocystdens) + assay +  (pf16c) + (pf25c) + (1|cluster),
  data = feed %>% mutate(assay = factor(assay, levels = c("MFA","DFA"))),
  family = "binomial")
prev_mod <- update(prev_mod,start=getME(prev_mod,c("theta","fixef")),control=glmerControl(optCtrl=list(maxfun=1e6)))

nms <- setNames(
  c("Monthly\nParasite\nPrevalence","Oocyst\nDensity","Direct\nFeeding\nAssay", "Pf16 Ct", "Pf25 Ct"),
  rownames(summary(prev_mod)$coefficients)[-1]
)
prev_or <- or_plot(prev_mod, 3.25, c(0.5,0.5), c(0.5,1,1.5,2, 2.5), nms, 1:5) +
  xlab("Adjusted Odds Ratio (aOR)")


# Bring plots all together
fig2a <- cowplot::plot_grid(pf16 + theme(legend.position = "none"),
                   pf25 + theme(legend.position = "none"),
                   cowplot::get_legend(int_prev),
                   ncol = 3, rel_widths = c(1,1,0.2),
                   labels = c("a", "b", "")) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

fig2b <- cowplot::plot_grid(int_prev + theme(legend.position = c(0.1,0.85), legend.background = element_rect(color = "black")),
                            prev_or,
                            ncol = 2, rel_widths = c(1,1.2),
                            labels = c("c", "d")) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

fig2 <- cowplot::plot_grid(fig2a, NA, fig2b, ncol = 1, rel_heights = c(1.2,0.05,1.2)) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

save_figs("oocyst_patterns", fig2, width = 12, height = 10)

### DEPRECATED

lm16 <- glm(pf16c ~ log(gfemia_mbf) + assay, data = feed %>% filter(gfemia_mbf>0))
lm25 <- glm(pf25c ~ log(gfemia_mbf) + assay, data = feed %>% filter(gfemia_mbf>0))
lmrev <- glm(log(gfemia_mbf) ~ pf16c + pf25c + assay + agecateg2, data = feed %>% filter(gfemia_mbf>0))



feed$gf_pred_16 <- 0
feed$gf_pred_16[feed$assay == "DFA"] <- (feed$pf16c[feed$assay == "DFA"] - lm16$coefficients[1])/lm16$coefficients[2]
feed$gf_pred_16[feed$assay == "MFA"] <- (feed$pf16c[feed$assay == "MFA"] - lm16$coefficients[1] - lm16$coefficients[3])/lm16$coefficients[2]

feed$gf_pred_25 <- 0
feed$gf_pred_25[feed$assay == "DFA"] <- (feed$pf25c[feed$assay == "DFA"] - lm25$coefficients[1])/lm25$coefficients[2]
feed$gf_pred_25[feed$assay == "MFA"] <- (feed$pf25c[feed$assay == "MFA"] - lm25$coefficients[1] - lm25$coefficients[3])/lm25$coefficients[2]


feed %>% filter(gf_pred_16>0) %>% ggplot(aes(pf16c, gf_pred_16)) + geom_point() + geom_smooth(method = "lm")
feed %>% filter(gf_pred_25<10) %>% ggplot(aes(exp(altpred), oocystdens, weight = gutdissected)) + geom_point() +
  geom_smooth(method="lm", aes(color="Exp Model"), formula= (y ~ exp(x)), se=FALSE, linetype = 1) +
  geom_smooth(span = 1)

feed$altpred <- predict(lmrev, type = "response", newdata = feed)
feed %>% filter(altpred<6) %>% ggplot(aes(exp(altpred), oocystdens)) + geom_point() +
  geom_smooth(method="lm", aes(color="Exp Model"), formula= (y ~ exp(x)), se=FALSE, linetype = 1) +geom_smooth(method = "lm")

feed %>% filter(gfemia_mbf>0) %>% ggplot(aes(log(gfemia_mbf), pf16c)) + geom_point() + geom_smooth(method = "lm")
feed %>% filter(gfemia_mbf>0) %>% ggplot(aes(log(gfemia_mbf), pf25c)) + geom_point() + geom_smooth(method = "lm")
feed %>% filter(altpred>0) %>% ggplot(aes(exp(altpred), oocystprevalence)) + geom_point() + geom_smooth(method = "lm") + xlim(c(0,1000)) + coor


feed %>% filter(altpred>0) %>% ggplot(aes(exp(altpred), oocystprevalence/100)) + geom_point() + binomial_smooth() + xlim(1,1000) +
  geom_smooth() +
  coord_trans(x="log10")


feed %>% filter(altpred>0) %>% ggplot(aes(exp(altpred), oocystprevalence/100)) + geom_point() + binomial_smooth() + xlim(c(0,500))

feed %>% filter(gf_pred_16>0) %>% ggplot(aes(exp(gf_pred_16), oocystprevalence/100)) + geom_point() + binomial_smooth() +
  scale_x_log10()

feed %>% filter(gf_pred_25>0) %>% ggplot(aes(exp(gf_pred_25), oocystdens)) + geom_point() + binomial_smooth() +
  scale_x_log10()
