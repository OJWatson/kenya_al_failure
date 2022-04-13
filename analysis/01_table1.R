# ------------------------------------------------------------------------------
# Pre. Read in our individual data
# ------------------------------------------------------------------------------

library(lme4)
all <- readRDS(cp_path("analysis/data/derived/all.rds"))
feed <- readRDS(cp_path("analysis/data/derived/feed.rds"))
all$cut_pfemia_mbf <- cut(all$pfemia_mbf, c(0,1,1000, 10000, 100000), include.lowest = TRUE)
all$thromb <- factor(c("No Thromb.", "Thromb.")[all$thromb+1], levels = c("No Thromb.", "Thromb."))
all$anemia <- factor(c("No Anemia", "Anemia")[all$anemia+1], levels = c("No Anemia", "Anemia"))

# ------------------------------------------------------------------------------
# 1. Table 1
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 1a. Create univariate analyses
# ------------------------------------------------------------------------------

asex_mod1 <- lme4::glmer(pf ~ 1 + gender + (1|stringid) + (1|homestead) + (1|cluster), all, family = "binomial")
asex_mod2 <- lme4::glmer(pf ~ 1 + agecateg2 + (1|stringid) + (1|homestead) + (1|cluster), all, family = "binomial")
asex_mod3 <- lme4::glmer(pf ~ 1 + cut_pfemia_mbf + (1|stringid) + (1|homestead) + (1|cluster), all, family = "binomial")
asex_mod4 <- lme4::glmer(pf ~ 1 + thromb + (1|stringid) + (1|homestead) + (1|cluster), all, family = "binomial")
asex_mod5 <- lme4::glmer(pf ~ 1 + anemia + (1|stringid) + (1|homestead) + (1|cluster), all, family = "binomial")
asex_mod6 <- lme4::glmer(pf ~ 1 + sickler + (1|stringid) + (1|homestead) + (1|cluster), all, family = "binomial")
asex_mods <- list(asex_mod1, asex_mod2, asex_mod3, asex_mod4, asex_mod5, asex_mod6)

sex_mod1 <- lme4::glmer(gametocyte ~ 1 + gender + (1|stringid) + (1|homestead) + (1|cluster), all, family = "binomial")
sex_mod2 <- lme4::glmer(gametocyte ~ 1 + agecateg2 + (1|stringid) + (1|homestead) + (1|cluster), all, family = "binomial")
sex_mod3 <- lme4::glmer(gametocyte ~ 1 + cut_pfemia_mbf + (1|stringid) + (1|homestead) + (1|cluster), all, family = "binomial")
sex_mod4 <- lme4::glmer(gametocyte ~ 1 + thromb + (1|stringid) + (1|homestead) + (1|cluster), all, family = "binomial")
sex_mod5 <- lme4::glmer(gametocyte ~ 1 + anemia + (1|stringid) + (1|homestead) + (1|cluster), all, family = "binomial")
sex_mod6 <- lme4::glmer(gametocyte ~ 1 + sickler + (1|stringid) + (1|homestead) + (1|cluster), all, family = "binomial")
sex_mods <- list(sex_mod1, sex_mod2, sex_mod3, sex_mod4, sex_mod5, sex_mod6)

# ------------------------------------------------------------------------------
# 1b. Summarise results into table
# ------------------------------------------------------------------------------

create_eff_ors <- function(mod, ref = "Female") {

  effs <- mod %>%
    broom.mixed::tidy(conf.int = TRUE) %>%
    mutate(across(c("estimate", "conf.low", "conf.high"), .fns = exp)) %>%
    mutate(across(c("estimate", "conf.low", "conf.high"), .fns = round, digits = 2)) %>%
    mutate(across(c("p.value"), .fns = round, digits = 3)) %>%
    select(term, estimate, conf.low, conf.high, p.value) %>%
    filter(!grepl("Intercept|sd", term)) %>%
    group_by(term) %>%
    summarise(or = paste0(estimate, " [", conf.low, "-", conf.high, "]"),
              p = p.value) %>%
    mutate(p = replace(p, which(p < 0.001), rep("<0.001", length(which(p < 0.001)))))

  effs2 <- rbind(data.frame(term = ref, or = 1, p = ""), effs)

  return(effs2)

}

refs <- c("Female", "5 years and under", "[0,1]", "No Thromb.", "No Anemia", "AA")

sex_effs <- do.call(rbind, lapply(seq_along(sex_mods), function(x) {create_eff_ors(sex_mods[[x]], refs[x])})) %>%
  mutate(term = gsub("gender|agecateg2|cut_pfemia_mbf|sickler|thromb|anemia", "", term))

asex_effs <- do.call(rbind, lapply(seq_along(asex_mods), function(x) {create_eff_ors(asex_mods[[x]], refs[x])})) %>%
  mutate(term = gsub("gender|agecateg2|cut_pfemia_mbf|sickler|thromb|anemia", "", term))

n1 <- rbind(
  all %>% group_by(gender) %>% summarise(n = n()) %>% setNames(c("term", "N")),
  all %>% group_by(agecateg2) %>% summarise(n = n()) %>% setNames(c("term", "N")),
  all %>% group_by(cut_pfemia_mbf) %>% summarise(n = n()) %>% setNames(c("term", "N")),
  all %>% group_by(thromb) %>% summarise(n = n()) %>% setNames(c("term", "N")),
  all %>% group_by(anemia) %>% summarise(n = n()) %>% setNames(c("term", "N")),
  all %>% group_by(sickler) %>% summarise(n = n()) %>% setNames(c("term", "N"))
) %>% na.omit()

# bring together and order the factors correctly
effs <- left_join(n1, asex_effs) %>% setNames(c("term", "N","sexor", "sexp")) %>%
  left_join(sex_effs) %>%
  setNames(c("Term", "N", "Asexual OR", "Asexual p-value","Sexual OR", "Sexual p-value"))

write.csv(effs, cp_path("analysis/tables/Table1.csv"))

##########


# ------------------------------------------------------------------------------
# 2. Table 2
# ------------------------------------------------------------------------------

mean_95ci <- function(x, n){
  x <- na.omit(x)
  n <- na.omit(n)
  bc <- round(Hmisc::binconf(sum(x), sum(n))*100, 1)
  paste0(bc[1], " [", bc[2], " - ", bc[3], "]")
}

median_iqr <- function(x){

  paste0(round(median(x, na.rm = TRUE),1), " [",
         paste0(round(quantile(x, na.rm = TRUE, probs = c(0.25,0.75)),1), collapse = " - "),
         "]")

}

# create summaries
tbl2 <- feed %>% group_by(assay) %>%
  summarise(n = n(),
            pmgd = mean_95ci(oocystposneg, as.integer(oocystposneg>=0)),
            ooc_prev = mean_95ci(posguts, gutdissected),
            ooc_dens = median_iqr(oocystdens),
            pf16ct = median_iqr(pf16c),
            pf25ct = median_iqr(pf25c)
            ) %>% t()

# p-values from simple t.test and wilcox due to random effects not being identifiable with reduced samples in feeding
tbl2ps <- feed %>%
  tableone::CreateTableOne(vars = c("oocystposneg", "oocystprevalence", "oocystdens", "pf16c", "pf25c"),
                           strata = "assay" , data = .) %>%
  print(nonnormal = c("oocystdens"), minMax = TRUE)

tbl2 <- cbind(tbl2, c("",tbl2ps[,3]))
write.csv(tbl2, cp_path("analysis/tables/Table2.csv"))

# ------------------------------------------------------------------------------
# 3. Table 3
# ------------------------------------------------------------------------------

# create summaries
tbl3a <- rbind(
  feed %>% select(subpf, oocystposneg) %>% na.omit %>%
  group_by(subpf) %>%
  summarise(is = mean_95ci(oocystposneg, as.integer(oocystposneg>=0))) %>%
    rename(term = subpf),
  feed %>% select(subgf, oocystposneg) %>% na.omit %>%
    group_by(subgf) %>%
    summarise(is = mean_95ci(oocystposneg, as.integer(oocystposneg>=0))) %>%
    rename(term = subgf),
  feed %>% select(pf_mono, oocystposneg) %>% na.omit %>%
    group_by(pf_mono) %>%
    summarise(is = mean_95ci(oocystposneg, as.integer(oocystposneg>=0))) %>%
    rename(term = pf_mono),
  feed %>% select(agecateg2, oocystposneg) %>% na.omit %>%
    group_by(agecateg2) %>%
    summarise(is = mean_95ci(oocystposneg, as.integer(oocystposneg>=0))) %>%
    rename(term = agecateg2)
  )

tbl3b <- rbind(
  feed %>% select(subpf, posguts, gutdissected) %>% na.omit %>%
    group_by(subpf) %>%
    summarise(is = mean_95ci(posguts, gutdissected)) %>%
    rename(term = subpf),
  feed %>% select(subgf, posguts, gutdissected) %>% na.omit %>%
    group_by(subgf) %>%
    summarise(is = mean_95ci(posguts, gutdissected)) %>%
    rename(term = subgf),
  feed %>% select(pf_mono, posguts, gutdissected) %>% na.omit %>%
    group_by(pf_mono) %>%
    summarise(is = mean_95ci(posguts, gutdissected)) %>%
    rename(term = pf_mono),
  feed %>% select(agecateg2, posguts, gutdissected) %>% na.omit %>%
    group_by(agecateg2) %>%
    summarise(is = mean_95ci(posguts, gutdissected)) %>%
    rename(term = agecateg2)
)

succ_mod1 <- glm(oocystposneg ~ 1 + subpf, feed, family = "binomial")
succ_mod2 <- glm(oocystposneg ~ 1 + subgf, feed, family = "binomial")
succ_mod3 <- glm(oocystposneg ~ 1 + pf_mono, feed, family = "binomial")
succ_mod4 <- glm(oocystposneg ~ 1 + agecateg2, feed, family = "binomial")
succ_mods <- list(succ_mod1, succ_mod2, succ_mod3, succ_mod4)

prev_mod1 <- glm(oocystprevalence/100 ~ 1 + subpf, feed, family = "binomial")
prev_mod2 <- glm(oocystprevalence/100 ~ 1 + subgf, feed, family = "binomial")
prev_mod3 <- glm(oocystprevalence/100 ~ 1 + pf_mono, feed, family = "binomial")
prev_mod4 <- glm(oocystprevalence/100 ~ 1 + agecateg2, feed, family = "binomial")
prev_mods <- list(prev_mod1, prev_mod2, prev_mod3, prev_mod4)

refs <- c("Microscopic asexual parasitaemia",
          "Microscopic gameotcytemia",
          "Pf + non-Pf infections",
          "<6 years")

succ_effs <- do.call(rbind, lapply(seq_along(succ_mods), function(x) {create_eff_ors(succ_mods[[x]], refs[x])})) %>%
  mutate(term = gsub("pf_mono|agecateg2", "", term)) %>%
  mutate(term = gsub("subpf","Submicroscopic asexual parasitemia", term)) %>%
  mutate(term = gsub("subgf","Submicroscopic gametocytemia", term)) %>%
  mutate(term = gsub("TRUE","Pf mono-infection", term))

prev_effs <- do.call(rbind, lapply(seq_along(prev_mods), function(x) {create_eff_ors(prev_mods[[x]], refs[x])})) %>%
  mutate(term = gsub("pf_mono|agecateg2", "", term)) %>%
  mutate(term = gsub("subpf","Submicroscopic asexual parasitemia", term)) %>%
  mutate(term = gsub("subgf","Submicroscopic gametocytemia", term)) %>%
  mutate(term = gsub("TRUE","Pf mono-infection", term))

# denominators
n1 <- rbind(
  feed %>% group_by(subpf) %>% summarise(n = n()) %>% setNames(c("term", "N")),
  feed %>% group_by(subgf) %>% summarise(n = n()) %>% setNames(c("term", "N")),
  feed %>% group_by(pf_mono) %>% summarise(n = n()) %>% setNames(c("term", "N")),
  feed %>% group_by(agecateg2) %>% summarise(n = n()) %>% setNames(c("term", "N"))
) %>% na.omit()

tbl3a$term[1:7] <- succ_effs$term[1:7]
tbl3b$term[1:7] <- succ_effs$term[1:7]
n1$term[1:7] <- succ_effs$term[1:7]

tbl3 <- left_join(n1,succ_effs) %>%
  left_join(tbl3a) %>%
  select(term, N, is, or, p) %>%
  setNames(c("term","N","succ_is","succ_or","succ_p")) %>%
  left_join(left_join(prev_effs, tbl3a) %>% select(term, is, or, p)) %>%
  setNames(c("Term", "N", "Infection Success","OR", "p-value","Infection Prevalence","OR", "p-value"))

write.csv(tbl3, cp_path("analysis/tables/Table3.csv"))

# ------------------------------------------------------------------------------
# 4. Supplementary Table 2
# ------------------------------------------------------------------------------

# write these results to supp table 2

a <- all %>%
  mutate(mbfpf = replace(mbfpf, pfpcr == 0, 0)) %>%
  mutate(mbfgam = replace(mbfgam, gfpcr == 0, 0)) %>%
  mutate(subpf = replace(subpf, is.na(subpf), FALSE)) %>%
  select(gender, age, agecateg2, mbfpf, pfpcr, subpf, mbfgam, gfpcr, subgf) %>%
  mutate(across(mbfpf:subgf, .fns = as.logical)) %>%
  tableone::CreateTableOne(vars = names(.)[-which(names(.) == "agecateg2")], strata = "agecateg2" , data = ., test = FALSE) %>%
  print(., formatOptions = list(big.mark = ","), nonnormal = "age", noSpaces = TRUE, minMax = TRUE)

b <- all %>%
  mutate(mbfpf = replace(mbfpf, pfpcr == 0, 0)) %>%
  mutate(mbfgam = replace(mbfgam, gfpcr == 0, 0)) %>%
  mutate(subpf = replace(subpf, is.na(subpf), FALSE)) %>%
  select(gender, age, agecateg2, mbfpf, pfpcr, subpf, mbfgam, gfpcr, subgf) %>%
  mutate(across(mbfpf:subgf, .fns = as.logical)) %>%
  tableone::CreateTableOne(vars = names(.)[-which(names(.) == "agecateg2")], data = ., test = FALSE) %>%
  print(., formatOptions = list(big.mark = ","), nonnormal = "age", noSpaces = TRUE, minMax = TRUE)

ab <- cbind(a, b)

# sens/spec

pf_sens_spec <- rbind(
  all %>% group_by(agecateg2) %>%
    summarise(sens_pf = caret::confusionMatrix(as.factor(mbfpf), as.factor(pfpcr), positive = "1")$byClass["Sensitivity"],
              spec_pf = caret::confusionMatrix(as.factor(mbfpf), as.factor(pfpcr), positive = "1")$byClass["Specificity"]),
  all %>% mutate(agecateg3 = "All") %>% group_by(agecateg3) %>%
    summarise(sens_pf = caret::confusionMatrix(as.factor(mbfpf), as.factor(pfpcr), positive = "1")$byClass["Sensitivity"],
              spec_pf = caret::confusionMatrix(as.factor(mbfpf), as.factor(pfpcr), positive = "1")$byClass["Specificity"]) %>%
    rename(agecateg2 = agecateg3)
) %>% t

gam_sens_spec <- rbind(
  all %>% group_by(agecateg2) %>%
    summarise(sens_gam = caret::confusionMatrix(as.factor(mbfgam), as.factor(gfpcr), positive = "1")$byClass["Sensitivity"],
              spec_gam = caret::confusionMatrix(as.factor(mbfgam), as.factor(gfpcr), positive = "1")$byClass["Specificity"]),
  all %>% mutate(agecateg3 = "All") %>% group_by(agecateg3) %>%
    summarise(sens_gam = caret::confusionMatrix(as.factor(mbfgam), as.factor(gfpcr), positive = "1")$byClass["Sensitivity"],
              spec_gam = caret::confusionMatrix(as.factor(mbfgam), as.factor(gfpcr), positive = "1")$byClass["Specificity"]) %>%
    rename(agecateg2 = agecateg3)
) %>% t


# bring all together

tab_supp_2 <- rbind(ab[1:6,], pf_sens_spec[-1,], ab[7:9,], gam_sens_spec[-1,])

write.csv(
  tab_supp_2,
  cp_path("analysis/tables/supp_table2.csv")
)

# ------------------------------------------------------------------------------
# 5. Supplementary Table 3
# ------------------------------------------------------------------------------

# create our grouped proportions
st3 <- rbind(
  all %>%
    mutate(gender = factor(gender, levels = c("Male", "Female"))) %>%
    group_by(gender) %>% summarise(n = n(), m = sum(mbfpf & pf), pm = m/n, q = sum(pf), pq = q/n) %>%
    rename(g = gender),
  all %>%
    group_by(agecateg2) %>% summarise(n = n(), m = sum(mbfpf & pf), pm = m/n, q = sum(pf), pq = q/n) %>%
    rename(g = agecateg2),
  all %>%
    summarise(n = n(), m = sum(mbfpf & pf), pm = m/n, q = sum(pf), pq = q/n) %>%
    mutate(g = "all", .before = 1)
)

# and check for significant difference from random effect models
p1 <- summary(lme4::glmer(mbfpf ~ 1 + gender + (1|stringid) + (1|homestead) + (1|cluster), all %>% mutate(mbfpf = replace(mbfpf, pfpcr == 0, 0)), family = "binomial"))$coefficients
p2 <- summary(lme4::glmer(mbfpf ~ 1 + agecateg2 + (1|stringid) + (1|homestead) + (1|cluster), all %>% mutate(mbfpf = replace(mbfpf, pfpcr == 0, 0)), family = "binomial"))$coefficients

p3 <- summary(lme4::glmer(pf ~ 1 + gender + (1|stringid) + (1|homestead) + (1|cluster), all %>% mutate(mbfpf = replace(mbfpf, pfpcr == 0, 0)), family = "binomial"))$coefficients
p4 <- summary(lme4::glmer(pf ~ 1 + agecateg2 + (1|stringid) + (1|homestead) + (1|cluster), all %>% mutate(mbfpf = replace(mbfpf, pfpcr == 0, 0)), family = "binomial"))$coefficients

pr <- function(p) {
  ifelse(min(tail(p[,4],-1)) < 0.001, "p<0.001", paste0("p = ", round(min(tail(p[,4],-1)),4)))
}

tab_supp_3 <- st3 %>%
  mutate(Group = paste0(g, " (N = ", n, ")"),
         Microscopy = paste0(m, " (", round(pm*100,1), "%)"),
         `p-value_Microscopy` = c(pr(p1), "", pr(p2),rep("",4)),
         qPCR = paste0(q, " (", round(pq*100,1), "%)"),
         `p-value_qPCR` = c(pr(p3), "", pr(p4),rep("",4))) %>%
  select(-g,-q,-n,-m,-pm,-pq)

write.csv(
  tab_supp_3,
  cp_path("analysis/tables/supp_table3.csv")
)

# ------------------------------------------------------------------------------
# 6. Supplementary Table 4s
# ------------------------------------------------------------------------------


# 4a sex numbers

all %>% select(agecateg2, mbfpf, pfpcr, pf16c, pf25c) %>% filter(pfpcr == 1) %>% mutate(marker = ifelse(mbfpf == 1 & pfpcr == 1, "M", "SM")) %>%
  group_by(agecateg2, marker) %>%
  summarise(n = n(),
            pf1625 = sum(!is.na(pf16c) & !is.na(pf25c), na.rm = TRUE),
            pf16 = sum(!is.na(pf16c) & is.na(pf25c), na.rm = TRUE),
            pf25 = sum(is.na(pf16c) & !is.na(pf25c), na.rm = TRUE),
            pf16o25 = sum(!is.na(pf16c) | !is.na(pf25c), na.rm = TRUE))

# 4a
tab_supp_4a <- rbind(all %>%
                       group_by(agecateg2) %>%
                       summarise(n = n(),
                                 pf1625 = sum(gametstage == "both", na.rm = TRUE),
                                 pf16 = sum(gametstage == "pf16", na.rm = TRUE),
                                 pf25 = sum(gametstage == "pf25", na.rm = TRUE),
                                 pf16opf25 = sum(pf1625, pf16, pf25, na.rm=TRUE)),
                     all %>%
                       summarise(n = n(),
                                 pf1625 = sum(gametstage == "both", na.rm = TRUE),
                                 pf16 = sum(gametstage == "pf16", na.rm = TRUE),
                                 pf25 = sum(gametstage == "pf25", na.rm = TRUE),
                                 pf16opf25 = sum(pf1625, pf16, pf25, na.rm=TRUE)) %>%
                       mutate(agecateg2 = "Total", .before =1)
) %>%
  mutate(across(pf1625:pf16opf25, .fns = ~ paste0(.x, " (",round(.x/n*100,1),"%)"))) %>%
  t()

write.csv(
  tab_supp_4a,
  cp_path("analysis/tables/supp_table4a.csv")
)

# 4b Sex numbers by asexual status

tab_supp_4b <- rbind(
  all %>% select(agecateg2, mbfpf, pfpcr, pf16c, pf25c, maturegamet,gametocyte, gfpcr, gametstage) %>%
  mutate(marker = ifelse(mbfpf == 1 & pfpcr == 1, "M", NA)) %>%
  mutate(marker = replace(marker, which(mbfpf == 0 & pfpcr == 1), "SM")) %>%
  filter(!is.na(marker)) %>%
  group_by(agecateg2, marker) %>%
  summarise(n = n(),
            pf1625 = sum(gametstage == "both", na.rm = TRUE),
            pf16 = sum(gametstage == "pf16", na.rm = TRUE),
            pf25 = sum(gametstage == "pf25", na.rm = TRUE),
            pf16opf25 = sum(gametocyte, na.rm=TRUE)),
  all %>% select(agecateg2, mbfpf, pfpcr, pf16c, pf25c, maturegamet,gametocyte, gfpcr, gametstage) %>%
    mutate(marker = ifelse(mbfpf == 1 & pfpcr == 1, "M", NA)) %>%
    mutate(marker = replace(marker, which(mbfpf == 0 & pfpcr == 1), "SM")) %>%
    filter(!is.na(marker)) %>%
    mutate(agecateg2 = "All") %>%
    group_by(agecateg2, marker) %>%
    summarise(n = n(),
              pf1625 = sum(gametstage == "both", na.rm = TRUE),
              pf16 = sum(gametstage == "pf16", na.rm = TRUE),
              pf25 = sum(gametstage == "pf25", na.rm = TRUE),
              pf16opf25 = sum(gametocyte, na.rm=TRUE))
  ) %>%
  mutate(pf16opf25 = paste0(pf16opf25, " (", round(pf16opf25/n*100, 1), "%)")) %>%
  mutate(n = paste0("N = ", n)) %>%
  t()

write.csv(
  tab_supp_4b,
  cp_path("analysis/tables/supp_table4b.csv")
)

# 4c sex numbers by age

tab_supp_4c <- rbind(
  all %>% filter(gametocyte == 1) %>% group_by(agecateg2) %>% summarise(m = mean(mbfpf)),
  all %>% filter(gametocyte == 1) %>% mutate(agecateg2 = "All") %>% group_by(agecateg2) %>% summarise(m = mean(mbfpf))
) %>%
  mutate(m = round(m*100,1)) %>%
  mutate(mn = 100-m) %>%
  t()

write.csv(
  tab_supp_4c,
  cp_path("analysis/tables/supp_table4c.csv")
)
