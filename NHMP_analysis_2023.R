#==================#
# NHMP Analysis 2023
#==================#
#
# System dependencies required for macOS users, install using Homebrew
# (https://brew.sh/):
#
#     brew install gdal udunits
#

# install.packages(c("devtools", "tidyverse"))
# devtools::install_github("inbo/camtraptor")
remove.packages("camtrapDensity")
devtools::install_github("MarcusRowcliffe/camtrapDensity")

library(tidyverse)
library(camtrapDensity)
source("helper_functions.R")

#===========================#
# DATAPACKAGE PROCESSING ####
#===========================#
## Read and merge datapackages ####
# pkgSp objects hold species-specific datapackages 
# as exported from agouti
# Sequence:
#   load data
#   check deployment models
#   merge list into a single datapackage
#   correct some inconsistencies in location names
#     - remove space before 23 (occurs in a few cases)
#     - remove number suffixes (all text after space - 
#       suffixes are inconsistent between species for some reason)
camtrapDir <- "./data/cameratrap/2023"
pkgHedgehog <- file.path(camtrapDir, "hedgehog") %>%
  list.dirs(recursive=FALSE) %>%
  file.path("datapackage.json") %>%
  map(read_camtrapDP) %>%
  map(check_deployment_models) %>%
  merge_camtrapDP() %>%
  mutate_camtrapDP(locationName = sub(" 23", "23", locationName),
                   locationName = sub(" .*", "", locationName),
                   obs = FALSE)
pkgBadger <- file.path(camtrapDir, "badger") %>%
  list.dirs(recursive=FALSE) %>%
  file.path("datapackage.json") %>%
  map(read_camtrapDP) %>%
  map(check_deployment_models) %>%
  merge_camtrapDP() %>%
  mutate_camtrapDP(locationName = sub(" 23", "23", locationName),
                   locationName = sub(" .*", "", locationName),
                   obs = FALSE)

# Merge species datapackages (add badger to hedgehog)
pkgAll <- pkgHedgehog %>%
  add_species(pkgBadger, "Meles meles") %>%
  mutate_camtrapDP(subsiteName = extract_subsite_names(locationName),
                   subsiteName = case_match(subsiteName,
                                            "[Redacted]" ~ "Durham",
                                            "[Redacted]" ~ "London 3",
                                            "[Redacted]" ~ "London 1",
                                            "[Redacted]" ~ "London 2",
                                            .default = subsiteName),
                   obs=FALSE)

## Deployment time corrections ####
# - Read deployment_dates spreadsheet
# - remove locations not present in datapackage
# - create reference date/times
# - add agouti deployment start date/time
deployment_corrections <- read.csv(file.path(camtrapDir, "deployment_dates.csv"),
                                             stringsAsFactors = FALSE) %>%
  filter(locationName %in% pkgAll$data$deployments$locationName) %>%
  mutate(true_start = dmy_hms(paste(true_start_date, true_start_time, tz)),
         true_end = dmy_hms(paste(true_end_date, true_end_time, tz))) %>%
  left_join(select(pkgAll$data$dep, locationName, start), by="locationName")
# Parse failure is OK (missing data, to be removed below)

# Remove location with no true deployment time data from datapackage and corrections table
missingLoc <- deployment_corrections %>%
  filter(is.na(true_start) | is.na(true_end)) %>%
  pull(locationName)
pkgAll <- subset_deployments(pkgAll, locationName != missingLoc)
deployment_corrections <- filter(deployment_corrections, locationName != missingLoc)

# Correct whole deployment times and update deployment start/end times
locMatch <- match(pkgAll$data$deployments$locationName, 
                  deployment_corrections$locationName)
pkgAll <- pkgAll %>%
  correct_time2(locName = deployment_corrections$locationName,
                wrongTime= deployment_corrections$start,
                rightTime = deployment_corrections$true_start) %>%
  mutate_camtrapDP(start = if_else(locationName %in% deployment_corrections$locationName,
                                   deployment_corrections$true_start[locMatch],
                                   start),
                   end = if_else(locationName %in% deployment_corrections$locationName,
                                 deployment_corrections$true_end[locMatch],
                                 end),
                   obs = FALSE)

# Check schedule (reordering by subsite)
pkgAll %>%
  mutate_camtrapDP(locationName = paste(subsiteName, substr(locationName, 1, 3)), obs=FALSE) %>%
  plot_deployment_schedule()

#============================#
# HABITAT DATA PROCESSING ####
#============================#
habDir <- "./data/habitat"
habitatAreas <- read.csv(file.path(habDir, "site_land_class_areas2024_long.csv")) %>%
  group_by(siteName, agg_className) %>%
  summarise(area = sum(area)) %>%
  ungroup()

## Define habitat clusters ####
# Elbow plot
clusData <- habitatAreas %>%
  pivot_wider(names_from = agg_className, values_from = area) %>%
  select(-siteName)
kk <- 2:10
twss <- sapply(kk, function(k) kmeans(clusData, k, nstart=10)$tot.withinss)
plot(kk, twss, type="b")

## Visualise habitat mix by cluster group ####
clust <- kmeans(clusData, 4, nstart=10)
habitatAreas %>%
  mutate(cluster = rep(clust$cluster, each=ncol(clust$centers))) %>%
  ggplot(aes(x = reorder(siteName, cluster),
             y = area, 
             fill=agg_className)) +
    geom_bar(position="fill", stat="identity") +
    annotate("text", 
             x = 1:length(unique(habitatAreas$siteName)), 
             y = 1.02,
             label = sort(clust$cluster)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## Map sites to dominant habitats ####
# Define cluster-specific dominant habitats 
habs <- habitatAreas %>%
  mutate(cluster = rep(clust$cluster, each=ncol(clust$centers))) %>%
  group_by(cluster, agg_className) %>%
  summarise(mnArea = mean(area)) %>%
  group_by(cluster) %>%
  filter(mnArea == max(mnArea)) %>%
  pull(agg_className) %>%
  case_match("Arable " ~ "Arable",
             "Broadleaf woodland " ~ "Woodland",
             "Built-up areas and gardens " ~ "Suburban",
             "Improved grassland " ~ "Pasture")
main_habitats <- data.frame(subsiteName = unique(habitatAreas$siteName),
                            habitat = habs[clust$cluster])
# Add habitats to datapackage
pkgAll <- pkgAll %>%
  mutate_camtrapDP(habitat = main_habitats$habitat[match(subsiteName,
                                                         main_habitats$subsiteName)])
main_habitats <- pkgAll$data$deployments %>%
  summarise(habitat = unique(habitat), .by=subsiteName)

#=================#
# SAVE POINT 1 ####
#=================#
save(pkgAll, file = file.path(camtrapDir, "pkgAll.Rdata"))
write.csv(main_habitats, file.path(habDir, "main_habitats.csv"),
          row.names = FALSE)

#=======================#
# DENSITY ESTIMATION ####
#=======================#
library(camtrapDensity)
library(tidyverse)
library(cowplot)
library(sbd)

## Read data ####
camtrapDir <- "./data/cameratrap/2023"
habDir <- "./data/habitat"
load(file.path(camtrapDir, "pkgAll.Rdata"))
main_habitats <- read.csv(file.path(habDir, "main_habitats.csv"))
          

## Summarise sample sizes ####
sp <- "Erinaceus europaeus"
sp <- "Meles meles"

ssize <- pkgAll$data$observations %>%
  left_join(select(pkgAll$data$deployments, deploymentID, subsiteName)) %>%
  summarise(n = sum(scientificName==sp, na.rm=TRUE), .by=subsiteName)
subsitesNonZero <- filter(ssize, n>0)$subsite
subsitesZero <- filter(ssize, n==0)$subsite

## fit effective detection radius models ####
rmod0 <- fit_detmodel(individualPositionRadius~1, pkgAll, sp, 
                      truncation = list(left=1, right=8))
plot(rmod0, pdf=T)
radii <- rmod0$edd %>% 
  as.data.frame() %>%
  mutate(lcl = estimate - 1.96*se,
         ucl = estimate + 1.96*se)


## fit effective detection angle models ####
amod0 <- fit_detmodel(individualPositionAngle~1, pkgAll, sp, unit="degree")
plot(amod0)
angles <- amod0$edd %>% 
  as.data.frame() %>%
  mutate(lcl = estimate - 1.96*se,
         ucl = estimate + 1.96*se)

## fit speed models ####
smod0 <- fit_speedmodel(pkgAll, sp, speed~1, pdf="l")
speeds <- smod0$estimate %>%
  rename(estimate = est)

## fit activity models ####
actmod0 <- fit_actmodel(pkgAll, species=sp, reps=1000)
actlevs <- actmod0@act %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  rename(estimate = act,
         lcl = 'lcl.2.5%',
         ucl = 'ucl.97.5%')

dranges <- data.frame(estimate = speeds$estimate * actlevs$estimate * 60^2*24/1000,
                      cv = sqrt(speeds$se/speeds$estimate)^2 +
                        (actlevs$se/actlevs$estimate)^2) %>%
  mutate(se = estimate * cv,
         lcl = exp(log(estimate) - 1.96*sqrt(log(1+cv^2))),
         ucl = exp(log(estimate) + 1.96*sqrt(log(1+cv^2))))

## fit density models ####
dmod0 <- rem_estimate(pkgAll, check_deployments=F, species=sp,
                      speed_model=smod0,
                      radius_model=rmod0,
                      angle_model=amod0,
                      activity_model=actmod0)
dmodS <- map(subsitesNonZero, \(s)
             subset_deployments(pkgAll, subsiteName == s) %>%
               rem_estimate(check_deployments=F, species=sp,
                            speed_model=smod0,
                            radius_model=rmod0,
                            angle_model=amod0,
                            activity_model=actmod0))

## Create density estimate table ####
dens <- map(c(list(dmod0), dmodS), \(m) 
            m$estimates %>%
              filter(rownames(.) == "density") %>%
              select(estimate, se, lcl95, ucl95)) %>%
  bind_rows() %>%
  mutate(subsiteName = c("OVERALL", subsitesNonZero)) %>%
  rename(lcl = lcl95,
         ucl = ucl95) %>%
  rbind(data.frame(estimate = 0, se=0, lcl = 0, ucl = 0, 
                   subsiteName = subsitesZero)) %>%
  left_join(ssize, by=join_by(subsiteName)) %>%
  left_join(main_habitats, by=join_by(subsiteName)) %>%
  mutate(n = ifelse(is.na(n), sum(n, na.rm=T), n)) %>%
  relocate(subsiteName, estimate)

#=================#
# SAVE POINT 2 ####
#=================#
write.csv(dens, paste0(sp,"_density.csv"), row.names = FALSE)

#==========#
# PLOTS ####
#==========#
library(ggimage)
library(tidyverse)
library(cowplot)


name_map <- c(
  "OVERALL" = "OVERALL",
  "[Redacted]" = "Durham 1",
  "[Redacted]" = "Lanarkshire",
  "[Redacted]" = "Lancashire",
  "[Redacted]" = "London 1",
  "[Redacted]" = "London 2",
  "[Redacted]" = "Nottinghamshire 1",
  "[Redacted]" = "Nottinghamshire 2",
  "[Redacted]" = "Dorset 1",
  "[Redacted]" = "Dorset 2",
  "[Redacted]" = "Isle of Wight",
  "[Redacted]" = "Hampshire",
  "[Redacted]" = "Leicestershire",
  "[Redacted]" = "London 3"
)


# Get data
densH <- read.csv(paste0("Erinaceus europaeus_density.csv"))
densB <- read.csv(paste0("Meles meles_density.csv"))

makePlot <- function(dens, label, img, 
                     zeroPos=0.01, zeroWidth=0.5,
                     imgSize=0.1, imgx=0.1, imgy=14){
  
  # define a *named* palette so mapping is stable
  base_cols <- c(
    Arable   = "#2078B3",
    Pasture  = "#6A3D9A",
    Woodland = "#B2DF8A",
    Suburban = "#FB9A99"
  )
  cols <- c(base_cols, OVERALL = "gray50")
  
  dens2 <- dens %>%
    dplyr::mutate(
      subsiteName = dplyr::recode(subsiteName, !!!name_map),
      habitat = factor(habitat, levels = names(base_cols)),   # keep legend order stable
      zero_all = (estimate == 0 & lcl == 0 & ucl == 0),
      
      estimate_plot = ifelse(estimate == 0, zeroPos, estimate),
      lcl_plot = ifelse(zero_all, NA_real_,
                        ifelse(lcl == 0, zeroPos / (1 + zeroWidth), lcl)),
      ucl_plot = ifelse(zero_all, NA_real_,
                        ifelse(ucl == 0, zeroPos * (1 + zeroWidth), ucl)),
      
      ord = ifelse(subsiteName == "OVERALL", 0,
                   rank(estimate_plot, ties.method = "random")),
      subsiteName = reorder(subsiteName, ord),
      
      # use a pseudo-habitat just for colouring the OVERALL row
      habitat_col = dplyr::if_else(subsiteName == "OVERALL", "OVERALL", as.character(habitat))
    )
  
  ggplot(dens2, aes(x = subsiteName, y = estimate_plot)) +
    geom_linerange(aes(ymin = lcl_plot, ymax = ucl_plot, colour = habitat_col),
                   lwd = 3, na.rm = TRUE) +
    geom_point(aes(colour = ifelse(estimate == 0, habitat_col, NA)),
               shape = 19, na.rm = TRUE) +
    # show only real habitats in the legend; map OVERALL via the named vector
    scale_colour_manual(
      name   = "Habitat",
      values = cols,                 # named mapping (including OVERALL)
      breaks = names(base_cols),     # legend shows only habitats
      na.value = "black"
    ) +
    geom_image(aes(image = img),
               data = tibble::tibble(subsiteName = imgy, estimate_plot = imgx),
               size = imgSize) +
    scale_y_continuous(breaks = c(zeroPos, 0.1, 1, 10, 100),
                       labels = c("0", "0.1", "1", "10", "100"),
                       transform = "log10") +
    labs(x = "", y = paste(label, "density")) +
    coord_flip() +
    theme_minimal()
}

plotH <- makePlot(densH, "Hedgehog", "Hedgehog.png")
plotB <- makePlot(densB, "Badger", "Badger.png", imgSize = 0.15)
legend <- get_legend(plotH)
plt <- plot_grid(plotlist = list(plotH + theme(legend.position = "none"),
                                 plotB + theme(legend.position = "none"),
                                 legend),
                 ncol = 3, rel_widths = c(2,2,1))
plt
ggsave("density plot.png", plt, height=5, width=7, bg="white")



