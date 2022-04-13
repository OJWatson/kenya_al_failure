# ------------------------------------------------------------------------------
# Pre. Read in our individual data and load packages
# ------------------------------------------------------------------------------

# Read in our data
library(tidyverse)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
library(rsatscan)
library(rgeoboundaries)
library(elevatr)
library(raster)

all <- readRDS(cp_path("analysis/data/derived/all.rds"))
geo <- all %>% select(homestead, lat, long, cluster) %>% unique

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
mycas2 <- mycas2 %>% filter(cluster %in% mygeo2$cluster) %>% as.data.frame()
mycont2 <- mycas2 %>% mutate(case  = as.integer(case != 1)) %>% filter(cluster %in% mygeo2$cluster)

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

# ------------------------------------------------------------------------------
# 2. Download relevant map shapes
# ------------------------------------------------------------------------------

# kenya shapes
kenya_bound <- rgeoboundaries::geoboundaries("Kenya")
kenya_bound2 <- rgeoboundaries::geoboundaries("Kenya", adm_lvl = "adm3") %>% filter(shapeName == "KISUMU RURAL")
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
            filter(P_VALUE < 0.05),
          aes(color = col), fill = NA) +
  scale_color_manual(name = "Hotspot Score\n", values = c("Coldspot" = "Blue", "Hotspot" = "Red")) +
  labs(x = "Longitude", y = "Latitude")

map <- cowplot::plot_grid(map1, map2, ncol = 1,
                          rel_heights = c(1,1),
                          align = "v", labels = c("a","b")) +
  theme(plot.background = element_rect(fill = "white"))
save_figs("map", map, 8, 10)
