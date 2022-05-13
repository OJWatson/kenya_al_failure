# ------------------------------------------------------------------------------
# Pre. Read in our individual data and load packages
# ------------------------------------------------------------------------------

# Read in our data
library(tidyverse)
library(rsatscan)
library(rgeoboundaries)
library(elevatr)
library(raster)

all <- readRDS(cp_path("analysis/data/derived/all.rds"))
geo <- all %>% dplyr::select(homestead, lat, long, cluster) %>% unique

# ------------------------------------------------------------------------------
# 1. SATSCAN
# ------------------------------------------------------------------------------

# aggregate data to homestead
mygeo2 <- geo[,-4] %>% group_by(homestead) %>%
  summarise(lat = median(lat), long = median(long)) %>%
  as.data.frame %>% mutate(homestead = as.integer(homestead)) %>%
  rename(cluster = homestead)
mycas2 <- all[,c("homestead", "pf", "visit")] %>%
  setNames(c("cluster", "case", "month")) %>%
  mutate(cluster = as.integer(cluster))
mycas2 <- mycas2 %>% dplyr::filter(cluster %in% mygeo2$cluster) %>% as.data.frame()
mycont2 <- mycas2 %>% mutate(case  = as.integer(case != 1)) %>% dplyr::filter(cluster %in% mygeo2$cluster)

# set up satscan runs
td <- file.path(tempdir(), "home")
dir.create(td)
write.geo(mygeo2, location = td, file = "mygeo")
write.cas(mycas2, location = td, file = "mycas")
write.ctl(mycont2, location = td, file = "mycont")
invisible(ss.options(reset=TRUE))
ss.options(list(ControlFile="mycont.ctl", CaseFile="mycas.cas"))
ss.options(list(StartDate="1", CoordinatesType=1, TimeAggregationUnits=4))
ss.options(list(EndDate="12", CoordinatesFile="mygeo.geo", ModelType=1, AnalysisType = 3, PrecisionCaseTimes=4))
ss.options(list(NonCompactnessPenalty=2, MaxTemporalSizeInterpretation=1, MaxTemporalSize=7, ScanAreas=3, UseDistanceFromCenterOption = "y",MaxSpatialSizeInDistanceFromCenter = 6.5))
ss.options(list(ProspectiveStartDate="1", ReportGiniClusters="n", LogRunToHistoryFile="n"))

write.ss.prm(td, "mybase")

# N.B. THIS WILL REQUIRE INSTALLING SATSCAN LOCALLY AND CHANGING THIS PATH
sat_res = satscan(td, "mybase",sslocation = "/home/oj/SaTScan")
supp_tab_5 <- sat_res$col %>%
  mutate(prevalence = OBSERVED/POPULATION) %>%
  dplyr::select(NUMBER_LOC, OBSERVED, POPULATION, RADIUS,
                prevalence, REL_RISK, P_VALUE) %>%
  mutate(P_VALUE = ifelse(P_VALUE < 0.0001, "<0.0001", round(P_VALUE, 4))) %>%
  mutate(RADIUS = round(RADIUS, 3)) %>%
  mutate(REL_RISK = round(REL_RISK, 2)) %>%
  mutate(prevalence = round(prevalence*100, 1)) %>%
  rename(Households = NUMBER_LOC,
         Cases = OBSERVED,
         Population = POPULATION,
         `Radius (km)` = RADIUS,
         Prevalence = prevalence,
         RR = REL_RISK,
         `p-value` = P_VALUE)

write.csv(supp_tab_5, cp_path("analysis/tables/supp_table5.csv"))

# ------------------------------------------------------------------------------
# 2. Download relevant map shapes
# ------------------------------------------------------------------------------

# kenya shapes
kenya_bound <- rgeoboundaries::geoboundaries("Kenya")
kenya_bound2 <- rgeoboundaries::geoboundaries("Kenya", adm_lvl = "adm3") %>% dplyr::filter(shapeName == "KISUMU RURAL")
kenya_bound3 <- rgeoboundaries::geoboundaries("Kenya", adm_lvl = "adm3")

# elevation data
elevation_data <- elevatr::get_elev_raster(locations = kenya_bound2, z = 11, clip = "locations")
elevation_data <- as.data.frame(elevation_data, xy = TRUE)
colnames(elevation_data)[3] <- "elevation"

# remove rows of data frame with one or more NA's,using complete.cases
elevation_data <- elevation_data[complete.cases(elevation_data), ]

# read in lakes
GLWD <- sf::read_sf(cp_path("analysis/data/raw/lake/GSHHS_c_L1.shp"))

# ------------------------------------------------------------------------------
# 3. Bring maps and prev and satscan together
# ------------------------------------------------------------------------------

# set up our various color palettes
col_palette <- colorspace::sequential_hcl(n = 7, h = c(36, 200), c = c(60, NA, 0), l = c(25, 95), power = c(0.7, 1.3))
col_palette <- c("#b4b4b4","#6d1708","#f4af01","#019617","#d4ea8c","#b1f0db")

# elevation map
map1 <- elevation_data %>%
  ggplot() +
  geom_sf(data = GLWD,
          mapping = aes(geometry = geometry),
          color = "grey",lwd=0.5,
          fill = "lightblue") +
  geom_sf(data = kenya_bound3, fill = "white", color = "black", lwd = 0.22) +
  geom_tile(aes(x = x, y = y, fill = elevation)) +
  geom_sf(data = kenya_bound2, color = "black", fill = NA, lwd = 0.2) +
  coord_sf() +
  scale_fill_gradientn(colors = col_palette, values = c(0,0.1,0.2,0.4, 0.6, 1), limits = c(1129, 1819), name = "Elevation (meters)\n") +
  scale_color_gradient2(midpoint=0.5, low="blue", mid="white",
                        high="red", space ="Lab" ) +
  xlim(c(34.4 , 34.7)) +
  ylim(c(-0.22,0.02)) +
  theme(axis.line = element_blank()) +
  labs(x = "Longitude", y = "Latitude")

# prevalence and hotspot map
pf_mean <- all %>% group_by(homestead, lat, long, cluster) %>%
  summarise(p = mean(pf))

map2 <- elevation_data %>%
  ggplot() +
  geom_sf(data = GLWD,
          mapping = aes(geometry = geometry),
          color = "grey",lwd=0.5,
          fill = "lightblue") +
  geom_sf(data = kenya_bound3, fill = "white", color = "black", lwd = 0.22) +
  geom_sf(data = kenya_bound2, color = "black", fill = NA) +
  geom_point(aes(long, lat, fill = p), shape = 21, stroke = 0.2, color = "black", data = pf_mean, inherit.aes = FALSE, size = 3, alpha  = 0.75) +
  coord_sf() +
  scale_fill_gradient2(midpoint=0.5, low="blue", mid="white",
                       high="red", space ="Lab", name = "Malaria Prevalence\n") +
  xlim(c(34.4 , 34.7)) +
  ylim(c(-0.21,0.02)) +
  geom_sf(data = sf::st_as_sf(sat_res$shapeclust) %>%
            mutate(col = c("Coldspot", "Hotspot")[as.integer(REL_RISK>1)+1]) %>%
            dplyr::filter(P_VALUE < 0.05),
          aes(color = col), fill = NA) +
  scale_color_manual(name = "Hotspot Score\n", values = c("Coldspot" = "Blue", "Hotspot" = "Red")) +
  labs(x = "Longitude", y = "Latitude")

map <- cowplot::plot_grid(map1, map2, ncol = 1,
                          rel_heights = c(1,1),
                          align = "v", labels = c("a","b")) +
  theme(plot.background = element_rect(fill = "white"))
save_figs("map", map, 8, 10)

# ------------------------------------------------------------------------------
# 4. Rainfall plot
# ------------------------------------------------------------------------------

rainfall <- read.csv(cp_path("analysis/data/raw/monthly_rainfall.csv")) %>%
  mutate(visit = lubridate::ceiling_date(as.Date(lubridate::dmonths(visit) + as.Date("2015-05-31")), "month"))

rainfall_gg <- all %>%
  group_by(visit) %>%
  summarise(Asexual = mean(pf),
            Gametocyte = mean(gametocyte)) %>%
  mutate(visit = lubridate::ceiling_date(as.Date(lubridate::dmonths(visit) + as.Date("2015-05-31")),"month")) %>%
  pivot_longer(Asexual:Gametocyte) %>%
  ggplot(aes(as.Date(visit), value, fill = name)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_smooth(aes(y = rainfall/max(rainfall), x = as.Date(visit), color = "Rainfall"),
              se = FALSE, span = 0.4,  data = rainfall, inherit.aes = FALSE) +
  ylim(c(0,1)) +
  scale_y_continuous(
    name = "Parasite Prevalence\n",
    sec.axis = sec_axis(~.*max(rainfall$rainfall), name="Monthly Rainfall (mm)\n", breaks = seq(0, 200, 20)),
    breaks = seq(0, 1, 0.1),
    labels = scales::label_percent()
  ) +
  theme_bw() +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10)) +
  scale_color_manual(values = "blue", name = "") +
  scale_fill_manual(values = viridis::cividis(2, begin = 0.4, end = 0.9), name = "") +
  scale_x_date(date_breaks = "months", date_labels = "%b-%y", limits = as.Date(c("2015-06-15","2016-06-15")), expand = c(0,0)) +
  xlab("")

save_figs("rainfall", rainfall_gg, 8, 4)
