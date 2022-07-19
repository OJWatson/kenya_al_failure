library(tidyverse)
devtools::load_all()
all <- readRDS(cp_path("analysis/data/derived/all.rds"))
feed <- readRDS(cp_path("analysis/data/derived/feed.rds"))

#-------------------------------------------------------------------------------
# Assign treatment Failures to the data
#-------------------------------------------------------------------------------

consec <- function(visit) {

  group <- rep(1, length(visit))
  for(i in seq_along(visit)[-1]) {

    if((visit[i]-visit[i-1]) != 1) {
      group[i] <- group[i-1]+1
    } else {
      group[i] <- group[i-1]
    }

  }

  return(group)
}

assign_tf <- function(pf) {

  tf <- rep(NA, length(pf))
  if(any(pf == 1)){
    # can't work out last in group so head -1
    inf <- which(head(pf,-1) == 1)
    for(i in seq_along(inf)) {
      if(pf[inf[i]+1] == 1) {
        tf[inf[i]] <- 1
      } else {
        tf[inf[i]] <- 0
      }
    }
  }

  tf[length(tf)] <- NA
  return(tf)
}

non_na_check <- function(x, b) {
  if(!is.na(x)) {
    x == b
  } else {
    return(FALSE)
  }
}
assign_transf <- function(pf_tf, oopn) {

  tf <- rep(NA, length(pf_tf))

  if(any(oopn == 1, na.rm = TRUE)){

    # can't work out first in group in group
    for(i in seq_along(oopn)[-1]) {
      if(!is.na(oopn[i]) & !is.na(pf_tf[i-1])){
        if(pf_tf[i-1] == 1) {
        if(oopn[i] == 1) {
          tf[i-1] <- 1
        } else {
          tf[i-1] <- 0
        }
        }
    }
    }
  }

  return(tf)
}

consecs <- all %>%
  filter(subjectid %in% unique(all %>% filter(pf == 1 | gametocyte == 1) %>% pull(subjectid))) %>%
  group_by(subjectid) %>%
  mutate(visit_grp = consec(visit)) %>%
  ungroup() %>%
  group_by(subjectid, visit_grp) %>%
  mutate(pf_tf = assign_tf(pf)) %>%
  mutate(gam_tf = assign_tf(gametocyte)) %>%
  mutate(pf_micro_tf = assign_tf(mbfpf)) %>%
  mutate(gam_micro_tf = assign_tf(mbfgam)) %>%
  mutate(transf = assign_transf(pf_tf, oocystposneg)) %>%
  mutate(cluster = substr(homestead, 1,2))

# proportion of feeds that were succesful and on individuals who failed the previous month
consecs$transf %>% table

#-------------------------------------------------------------------------------
# Model Treatment Failures
#-------------------------------------------------------------------------------

mod <- lme4::glmer(pf_tf ~ (age) + (visit_prev) + (1|homestead),
                   data = consecs %>% filter(pf == 1 & !is.na(pf_tf)) %>% mutate(visit_prev = visit_prev/100),
                   family = "binomial")
mod <- update(mod,start=lme4::getME(mod,c("theta","fixef")), control=lme4::glmerControl(optCtrl=list(maxfun=1e5)))

# explore the random effects
refp <- sjPlot::plot_model(mod, type="re", sort.est = TRUE)
refp_data <- refp$data %>%
  mutate(psig = ifelse(conf.low > 1, "High", ifelse(conf.high < 1, "Low", "None")))
refp_gg <- refp_data %>%
  ggplot() +
  geom_hline(yintercept = 1, color = "black") +
  geom_point(aes(y = estimate, x = term, color = psig)) +
  geom_errorbar(aes(ymax = conf.high, ymin = conf.low, x = term, color = psig)) +
  scale_y_log10(limits = c(0.1,10)) +
  scale_x_discrete(breaks = refp_data$term[refp_data$psig != "None"]) +
  scale_color_manual(values = c("High" = rb[1], "Low" = rb[2], "None" = "grey"), name = "Significant") +
  xlab("Homestead") +
  ylab("Homestead-level Random Effect for Treatment Failure") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10), panel.grid.minor.x = element_line(),
        legend.text = element_text(size = 10)) +
  coord_flip()

#-------------------------------------------------------------------------------
# Show the main drivers of treatment failure
#-------------------------------------------------------------------------------

pr1a <- ggeffects::ggpredict(mod, c("visit_prev [all]"))
pr1b <- plot(ggeffects::ggpredict(mod, c("visit_prev [all]", "homestead"), type = "random"))

pr1a_gg <- pr1a %>%
  mutate(psig = "Median") %>%
  ggplot(aes(x, predicted)) +
  # geom_line() +
  geom_line(aes(group = group, color = psig), data = left_join(pr1b$data, refp_data %>% mutate(group = term), by = "group"), alpha = 0.4) +
  geomtextpath::geom_textline(label = "Population-level", color = "black", lwd = 4) +
  scale_color_manual(
    values = c("High" = rb[1], "Low" = rb[2], "None" = "grey", "Overall" = "white"),
    labels = c("High" = "Significant Increased Risk (p<0.05)",
               "Low" = "Significant Decreased Risk (p<0.05)",
               "None" = "Non-significant Change in Risk (p>0.05)", "Overall" = ""),
    name = "Homestead-level Random Effects:") +
  xlab("Cluster-level Monthly Malaria Prevalence ") +
  ylab("Model Predicted Probability \nof Treatment Failure") +
  theme_bw() +
  theme(legend.position = "top") +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

pr2a <- ggeffects::ggpredict(mod, c("age [all]"))
pr2b <- plot(ggeffects::ggpredict(mod, c("age [all]", "homestead"), type = "random"))

pr2a_gg <- pr2a %>%
  mutate(psig = "Median") %>%
  ggplot(aes(x, predicted)) +
  # geom_line() +
  geom_line(aes(group = group, color = psig), data = left_join(pr2b$data, refp_data %>% mutate(group = term), by = "group"), alpha = 0.4) +
  geomtextpath::geom_textline(label = "Population-level", color = "black", lwd = 4) +
  scale_color_manual(
    values = c("High" = rb[1], "Low" = rb[2], "None" = "grey", "Overall" = "white"),
    labels = c("High" = "Significant Increased Risk (p<0.05)",
               "Low" = "Significant Decreased Risk (p<0.05)",
               "None" = "Non-significant Change in Risk (p>0.05)", "Overall" = ""),
    name = "Homestead-level Random Effects:") +
  xlab("Age (Years) ") +
  ylab("Model Predicted Probability \nof Treatment Failure") +
  theme_bw() +
  theme(legend.position = "top") +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))


#-------------------------------------------------------------------------------
# Create figure for dynamics of interesting groups
#-------------------------------------------------------------------------------

labeller <- as.list(consecs %>% ungroup %>%
                      filter(homestead %in% refp_data$term[refp_data$psig == "High"]) %>%
                      select(stringid, age) %>%
                      unique %>%  pull(age) %>%
                      paste0("L", round(.,0), " years") %>%
                      gsub("(.*L)(.*)","\\2",.))
names(labeller) <- consecs %>% ungroup %>%
  filter(homestead %in% refp_data$term[refp_data$psig == "High"]) %>%
  select(stringid, age) %>%
  unique %>%  pull(stringid)

for(i in seq_along(labeller)) {
  labeller[[i]] <- paste0(
    "Homestead: ", substr(names(labeller)[i],1,4),
    ", ID ", substr(names(labeller)[i],5,6),
    ": " ,labeller[[i]])
}

to_add <- data.frame("stringid" = c("070901",  "230305"), "age" = 99)
labeller <- c(labeller, setNames(as.list(paste0(to_add$stringid,": ", to_add$age)), to_add$stringid))

labeller_func <- function(variable, value){
  return(labeller[value])
}


# if they were chosen for feeding and gave rise to infections they must be gametoctye positive despite failing pcr/micro
# which we want to show in the plots but not in previous stats as technically not gam. positive by pcr
consecs <- consecs %>% mutate(gametocyte = replace(gametocyte, gametocyte == 0 & oocystposneg == 1, 1))
consecs <- consecs %>% mutate(gfpcr = replace(gfpcr, gfpcr == 0 & ((!is.na(pf16c) | !is.na(pf25c)) == 1), 1))

consecs$gam_stat <- ifelse(consecs$gametocyte == 1, "Gametocyte +ve, No Feeding Assay", "Gametocyte -ve")
consecs <- consecs %>%
  mutate(gam_stat = replace(gam_stat, gametocyte == 1 & !is.na(assay) & oocystposneg, "Gametocyte +ve, Mosquito Infected")) %>%
  mutate(gam_stat = replace(gam_stat, gametocyte == 1 & !is.na(assay) & !oocystposneg, "Gametocyte +ve, Mosquito Not Infected"))

fill_cols <- c("#0f7ba2", "#43b284", "#fab255", "#dd5129")
tf_patterns_gg <- consecs %>%
  filter(homestead %in% refp_data$term[refp_data$psig == "High"]) %>%
  full_join(data.frame("stringid" = c("070901", "230305"), "age" = 99)) %>%
  ggplot(aes(visit, plu, group = as.factor(visit_grp))) +
  geom_line(lwd = 0.25) +
  geom_point(aes(fill = as.factor(gam_stat), shape = as.factor(mbfpf)), stroke = 0.2, size = 3, na.rm = TRUE) +
  facet_wrap(~stringid, ncol = 5, labeller = labeller_func) +
  xlab("Visit Month") +
  ylab("Asexual Stage Ct Value") +
  theme_bw() +
  #scale_shape_manual(values = c(21:25), name = "Gametocyte Status") +
  scale_fill_manual(
    values = c("Gametocyte -ve" = fill_cols[1],
               "Gametocyte +ve, No Feeding Assay"=fill_cols[2],
               "Gametocyte +ve, Mosquito Not Infected" = fill_cols[3],
               "Gametocyte +ve, Mosquito Infected" = fill_cols[4]), name = "Gametocyte Status and \nFeeding Assay Outcome") +
  #scale_shape_manual(values = c("0"=21, "1" = 22), labels = c("1"="TRUE", "0"="FALSE"), name = "Gametocyte +ve") +
  scale_shape_manual(values = c("1"=21,"0"=22),
                    labels = c("1"="TRUE", "0"="FALSE"), name = "Asexual Microscopy +ve") +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  guides(fill=guide_legend(override.aes=list(shape=21, size = 4))) +
  guides(shape=guide_legend(override.aes=list(fill="white", stroke = 1, size = 4))) +
  theme(strip.background = element_rect(fill = "white"))
tf_patterns_gg

tf_gg <- cowplot::plot_grid(
  cowplot::plot_grid(
    cowplot::get_legend(pr1a_gg + theme(legend.text = element_text(size = 11),
                                        legend.title = element_text(size = 12))),
    cowplot::plot_grid(refp_gg + theme(legend.position = "none"),
                       pr1a_gg + theme(legend.position = "none"),
                       pr2a_gg + theme(legend.position = "none"),
                       ncol = 3, rel_widths = c(1,1,1),labels = "auto"),
    ncol = 1, rel_heights = c(0.1,1)),
  NA,
  tf_patterns_gg + theme(legend.position = "none"),
  ncol = 1, labels = c("","","d"),
  rel_heights = c(1,0.02, 1)
) + theme(plot.background = element_rect(fill = "white"))

tf_gg <- cowplot::plot_grid(tf_gg, cowplot::get_legend(tf_patterns_gg), rel_widths = c(1,0.2)) +
  theme(plot.background = element_rect(fill = "white"))
save_figs("treatment_failures", tf_gg, width = 16, height = 10)


#-------------------------------------------------------------------------------
# Create figure for dynamics for all groups
#-------------------------------------------------------------------------------

alldf <- consecs %>% ungroup %>%
  filter(visit_grp %in% names(which(table(visit_grp)>1))) %>%
  filter(stringid %in% (consecs %>% filter(pf == 1) %>% group_by(stringid) %>% summarise(n = n()) %>% filter(n > 0) %>% pull(stringid)))

labeller <- as.list(
  alldf %>%
    select(stringid, age) %>%
    unique %>%  pull(age) %>%
    paste0("L", round(.,0), " y") %>%
    gsub("(.*L)(.*)","\\2",.))

names(labeller) <- alldf %>%
  select(stringid, age) %>%
  unique %>%  pull(stringid)

for(i in seq_along(labeller)) {
  labeller[[i]] <- paste0(
    names(labeller)[i],
    ": " ,
    labeller[[i]])
}


labeller_func <- function(variable, value){
  return(labeller[value])
}

for(i in 1:5){

tf_patterns_gg <- alldf %>%
  ggplot(aes(visit, plu, group = as.factor(visit_grp))) +
  scale_fill_discrete() +
  geom_line(lwd = 0.25) +
  geom_point(aes(fill = as.factor(gam_stat), shape = as.factor(mbfpf)), stroke = 0.2, size = 2) +
  # facet_wrap(~stringid, ncol = 10, labeller = labeller_func) +
  ggforce::facet_wrap_paginate(~stringid, ncol = 7, nrow = 12, labeller = labeller_func, page = i) +
  xlab("Visit Month") +
  ylab("Asexual Stage Ct Value") +
  theme_bw() +
  #scale_shape_manual(values = c(21:25), name = "Gametocyte Status") +
  scale_fill_manual(
    values = c("Gametocyte -ve" = fill_cols[1],
               "Gametocyte +ve, No Feeding Assay"=fill_cols[2],
               "Gametocyte +ve, Mosquito Not Infected" = fill_cols[3],
               "Gametocyte +ve, Mosquito Infected" = fill_cols[4]), name = "Gametocyte Status and\nFeeding Assay Outcome:") +
  #scale_shape_manual(values = c("0"=21, "1" = 22), labels = c("1"="TRUE", "0"="FALSE"), name = "Gametocyte +ve") +
  scale_shape_manual(values = c("1"=21,"0"=22),
                     labels = c("1"="TRUE", "0"="FALSE"), name = "Asexual Microscopy +ve:") +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  guides(fill=guide_legend(override.aes=list(shape=21, size = 4), nrow = 2)) +
  guides(shape=guide_legend(override.aes=list(fill="white", stroke = 1, size = 4), nrow = 2)) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8),
        legend.position = "top")

save_figs(paste0("all_treatment_failures_",i), tf_patterns_gg, width = 10, height = 10)

}
qpdf::pdf_combine(
  list.files(cp_path("analysis/plots/"),
             pattern = "all_treatment_failures_\\d\\.pdf",
             full.names = TRUE),
  cp_path("analysis/plots/all_treatment_failures.pdf")
  )


# ------------
# Numbers
# ------------

# total infected at least once (asexual or sexual)
all %>% filter(stringid %in% (consecs %>% filter(gametocyte == 1 | pf == 1) %>% pull(stringid) %>% unique)) %>%
  pull(stringid) %>% unique %>% length

# two or more months sampling and gam
c(all %>%
  filter(stringid %in% (consecs %>% filter(gametocyte == 1) %>% pull(stringid) %>% unique)) %>%
  group_by(subjectid) %>%
  mutate(visit_grp = consec(visit)) %>%
  summarise(m = max(table(visit_grp))>1) %>%
  pull(m) %>% sum,
  all$subjectid %>% unique %>% length)

# two or more months sampling and pf
c(all %>%
    filter(stringid %in% (consecs %>% filter(pf == 1) %>% pull(stringid) %>% unique)) %>%
    group_by(subjectid) %>%
    mutate(visit_grp = consec(visit)) %>%
    summarise(m = max(table(visit_grp))>1) %>%
    pull(m) %>% sum,
  all$subjectid %>% unique %>% length)

# Tables of LPF episodes
consecs %>%
  filter(stringid %in% (consecs %>% filter(gametocyte == 1) %>% pull(stringid) %>% unique)) %>%
  select(subjectid, visit, visit_grp, pf, gametocyte, gam_tf) %>%
  ungroup %>%
  group_by(subjectid) %>%
  complete(visit = 1:12) %>%
  mutate(gam_tf = replace_na(gam_tf, 0)) %>%
  summarise(m = sum(diff(gam_tf)==-1)) %>%
  pull(m) %>% table

consecs %>%
  filter(stringid %in% (consecs %>% filter(pf == 1) %>% pull(stringid) %>% unique)) %>%
  select(subjectid, visit, visit_grp, pf, gametocyte, pf_tf) %>%
  filter(visit_grp %in% names(which(table(visit_grp)>1))) %>%
  ungroup %>%
  group_by(subjectid) %>%
  complete(visit = 1:12) %>%
  mutate(pf_tf = replace_na(pf_tf, 0)) %>%
  summarise(m = sum(diff(pf_tf)==-1)) %>%
  pull(m) %>% table

# Tables of LPF episodes

lpf_lengths <- function(tf) {

  if(!any(tf == 1, na.rm = TRUE)) {
    return(0)
  }

  tf[tf == 0] <- NA
  rl <- rle(is.na(tf))
  i1 <- rl$lengths>=1 & rl$values

  lst <- split(tf, rep(cumsum(c(TRUE, i1[-length(i1)])), rl$lengths)) %>%
    lapply(na.omit) %>%
    lengths() %>%
    max()
}

prepf <- consecs %>%
  filter(stringid %in% (consecs %>% filter(pf == 1) %>% pull(stringid) %>% unique)) %>%
  select(subjectid, visit, visit_grp, pf, gametocyte, pf_tf) %>%
  ungroup %>%
  group_by(subjectid) %>%
  complete(visit = 1:12)

pregf <- consecs %>%
  filter(stringid %in% (consecs %>% filter(gametocyte == 1) %>% pull(stringid) %>% unique)) %>%
  select(subjectid, visit, visit_grp, pf, gametocyte, gam_tf) %>%
  ungroup %>%
  group_by(subjectid) %>%
  complete(visit = 1:12)


lpfgg <- rbind(
  prepf %>% summarise(m = lpf_lengths(pf_tf)+1) %>% mutate(type = "Asexual Infection"),
  pregf %>% summarise(m = lpf_lengths(gam_tf)+1) %>% mutate(type = "Sexual Infection")) %>%
  ggplot(aes(x = m)) + geom_bar(stat = "count", fill = "white", color = "black") +
  scale_x_continuous(breaks = 1:12) +
  theme_bw() +
  facet_wrap(~type) +
  xlab("Duration in months of longest observed infection") +
  ylab("Count")

save_figs("lpf_duration", lpfgg, width = 8, height = 4)


# numbers for the fig:
cumsum((pregf %>% summarise(m = lpf_lengths(gam_tf)+1) %>% pull(m) %>% table))/
  sum((pregf %>% summarise(m = lpf_lengths(gam_tf)+1) %>% pull(m) %>% table))

cumsum((prepf %>% summarise(m = lpf_lengths(pf_tf)+1) %>% pull(m) %>% table))/
  sum((prepf %>% summarise(m = lpf_lengths(pf_tf)+1) %>% pull(m) %>% table))
