library(tidyverse)
library(haven)
library(rdhs)
devtools::load_all()

# ------------------------------------------------------------------------------
# 1. Read in our feed data
# ------------------------------------------------------------------------------

# Read in our raw feed data
feed <- haven::read_dta(cp_path("analysis/data/raw/903_feeding.dta"))

# date formats
feed <- feed %>%
  mutate(sampledate = as.Date(sampledate, "%d/%b/%y")) %>%
  mutate(feedingassaydate = as.Date(feedingassaydate, "%d-%b-%y")) %>%
  mutate(gutdissectiondate = as.Date(gutdissectiondate, "%d-%b-%y")) %>%
  mutate(glanddissectiondate = as.Date(glanddissectiondate, "%d-%b-%y")) %>%
  mutate(across(ends_with("c"),  ~ as.numeric(.x))) %>%
  mutate(age = floor((as.Date("2015-07-01") - dateofbirth)/365)) %>%
  select(-c( pf25u, pf16u, agecateg, latitude, longitude))

feed <- rdhs:::factor_format(feed, reformat = TRUE)$dataset %>%
  mutate(stringid = replace(stringid, nchar(stringid) == 5, paste0("0", stringid[nchar(stringid) == 5]))) %>%
  mutate(homestead = replace(homestead, nchar(homestead) == 3, paste0("0", homestead[nchar(homestead) == 3]))) %>%
  mutate(cluster = substr(homestead, 1, 2))

feed <- feed %>%
  mutate(uid = paste0(stringid, "_", visit)) %>%
  select(-(doorcode:sickler))

# ------------------------------------------------------------------------------
# 2. Read in our individual data and format
# ------------------------------------------------------------------------------

all <- haven::read_dta(cp_path("analysis/data/raw/all.dta"))
age_labels <- names(attr(all$agecateg2, "labels"))

all <- all %>%
  mutate(sampledate = datesample) %>%
  mutate(age = as.numeric((as.Date("2015-07-01") - dateofbirth)/365)) %>%
  mutate(slidedate = as.Date(slidedate, "%d-%b-%y")) %>%
  select(-c(visitdate, gender_c, datesample))

all <- rdhs:::factor_format(all, reformat = TRUE)$dataset %>%
  mutate(stringid = replace(stringid, nchar(stringid) == 5, paste0("0", stringid[nchar(stringid) == 5]))) %>%
  mutate(homestead = replace(homestead, nchar(homestead) == 3, paste0("0", homestead[nchar(homestead) == 3]))) %>%
  mutate(cluster = substr(homestead, 1, 2)) %>%
  mutate(agecut = cut(age, c(0,6,12,25,1000))) %>%
  mutate(agecut = factor(agecut, levels = age_labels)) %>%
  mutate(agecateg2 = factor(agecateg2, levels = age_labels))

all <- all %>%
  mutate(uid = paste0(stringid, "_", visit)) %>%
  mutate(across(c("pf", "gametocyte", "pfpcr", "gfpcr"), .fns = ~replace(.x, .x == "Positive", "1"))) %>%
  mutate(across(c("pf", "gametocyte", "pfpcr", "gfpcr"), .fns = ~replace(.x, .x == "Negative", "0"))) %>%
  mutate(across(c("pf", "gametocyte", "pfpcr", "gfpcr"), .fns = as.numeric)) %>%
  mutate(thromb = as.integer(plt < 120 )) %>%
  filter(!is.na(pf) & !is.na(gametocyte))

# visit dates wrong
all$visit[all$stringid == "410201" & all$visit == 6] <- 7

# ------------------------------------------------------------------------------
# 3. Merge across to bring over outcomes and fill in gaps
# ------------------------------------------------------------------------------

# merge in the geo data
geo <- haven::read_dta(cp_path("analysis/data/raw/homesteadgeo.dta")) %>%
  rename(lat = latitude, long = longitude) %>% group_by(homestead) %>%
  summarise(lat = median(lat), long = median(long)) %>%
  mutate(cluster = substr(homestead, 1, 2))

misshh <- data.frame(
  "homestead" = unique(all$homestead[which(!all$homestead %in% geo$homestead)])
) %>%
  mutate(cluster = substr(homestead, 1, 2))

geo <- rbind(geo,
             left_join(misshh,
                       geo %>% filter(cluster %in% misshh$cluster) %>%
                         group_by(cluster) %>%
                         summarise(lat = mean(lat),
                                   long = mean(long))
             ))

all <- all %>% group_by(visit, cluster) %>% mutate(visit_prev = mean(pf)*100) %>% ungroup

feed2 <- left_join(feed, all, by = "uid") %>%
  select(-ends_with(".x")) %>%
  rename_with(.fn = function(x){gsub(".y", "", x)}, .cols = ends_with(".y")) %>%
  mutate(pf_mono = pf == 1 & pm == "no" & po == "no") %>%
  left_join(na.omit(unique(geo[,c("lat","long","homestead")])), by = "homestead")

all2 <- left_join(all, na.omit(unique(geo[,c("lat","long","homestead")])), by = "homestead") %>%
  left_join(feed, by = "uid") %>%
  select(-ends_with(".y")) %>%
  rename_with(.fn = function(x){gsub(".x", "", x)}, .cols = ends_with(".x"))

# fill in blanks due to poor link where possible from the complete sampleid
feed2$subjectid[is.na(feed2$subjectid)] <- as.integer(gsub("-","",feed2$sampleid[is.na(feed2$subjectid)]))
feed2$stringid[is.na(feed2$stringid)] <- gsub("-","",feed2$sampleid[is.na(feed2$stringid)])
feed2$homestead[is.na(feed2$homestead)] <- substr(gsub("-","",feed2$sampleid[is.na(feed2$homestead)]),1,4)

# and fill in per person info
fill_missing <- function(feed2, all, var = "age") {
  feed2[[var]][is.na(feed2[[var]])] <- all[[var]][match(feed2$subjectid[is.na(feed2[[var]])], all$subjectid)]
  return(feed2)
}

to_fill <- c("doorcode", "gender", "maritalstatus","dateofbirth",
             "age", "agecateg", "sickle_type", "agecateg2",
             "sickler", "cluster", "lat", "long", "visit_prev")

for(i in to_fill) {
  feed2 <- fill_missing(feed2, all2, i)
}

feed2 <- feed2 %>% mutate(age = (as.Date(feedingassaydate) - as.Date(dateofbirth))/365) %>%
  mutate(assay = replace(assay, age < 7, "MFA"))

saveRDS(all2, cp_path("analysis/data/derived/all.rds"))
saveRDS(feed2, cp_path("analysis/data/derived/feed.rds"))

