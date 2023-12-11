library(terra)
library(SSDM)
library(readr)
library(dplyr)
library(usdm)
library(ggplot2)
library(rasterVis)
library(tidyr)
library(CoordinateCleaner)
library(spThin)

setwd("/home/borshon/infinity/jp_geo")

set.seed(13)

bioclims <- terra::rast("jp_R/bioclims_jp.tif")
f_245 <- terra::rast("jp_R/f_245_100jp.tif")
f_585 <- terra::rast("jp_R/f_585_100jp.tif")
jp <- terra::vect("JPN_adm/JPN_adm0.shp")
jp_w <- terra::vect("jp_water_3.shp")
carp_data <- read_csv("carp.csv")

carp_data <- select(carp_data, species = "species", long = "decimalLongitude",
                  lat = "decimalLatitude")


carp_d_cleaned <- cc_zero(carp_data, lon = "long", lat = "lat") %>%
  cc_dupl(lon = "long", lat = "lat", species = "species") %>%
  cc_sea(lon = "long", lat = "lat")


carp_d_thinned <- thin(loc.data = carp_d_cleaned, lat.col = "lat",
                       long.col = "long", spec.col = "species",
                       thin.par = 5, reps = 1, out.dir = getwd())


write_csv(carp_d_thinned, "jp_R/carp_thined.csv")

carp_2_thined <- read_csv("jp_R/carp_thined.csv")


carp_v <- vect(carp_2_thined, geom = c("long", "lat"), crs = crs(bioclims))

terra::plot(jp)
terra::points(carp_d_v, cex = 0.6, col="red", main = "Observation from gbif")


nms <- c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8",
         "bio9","bio10","bio11","bio12","bio13","bio14","bio15",
         "bio16","bio17","bio18","bio19")


jp <- terra::vect("JPN_adm/JPN_adm0.shp")

curr <- terra::rast(c("../worldclim/cur_bio/wc2.1_2.5m_bio_1.tif","../worldclim/cur_bio/wc2.1_2.5m_bio_2.tif",
                      "../worldclim/cur_bio/wc2.1_2.5m_bio_3.tif","../worldclim/cur_bio/wc2.1_2.5m_bio_4.tif",
                      "../worldclim/cur_bio/wc2.1_2.5m_bio_5.tif","../worldclim/cur_bio/wc2.1_2.5m_bio_6.tif",
                      "../worldclim/cur_bio/wc2.1_2.5m_bio_7.tif","../worldclim/cur_bio/wc2.1_2.5m_bio_8.tif",
                      "../worldclim/cur_bio/wc2.1_2.5m_bio_9.tif","../worldclim/cur_bio/wc2.1_2.5m_bio_10.tif",
                      "../worldclim/cur_bio/wc2.1_2.5m_bio_11.tif","../worldclim/cur_bio/wc2.1_2.5m_bio_12.tif",
                      "../worldclim/cur_bio/wc2.1_2.5m_bio_13.tif","../worldclim/cur_bio/wc2.1_2.5m_bio_14.tif",
                      "../worldclim/cur_bio/wc2.1_2.5m_bio_15.tif","../worldclim/cur_bio/wc2.1_2.5m_bio_16.tif",
                      "../worldclim/cur_bio/wc2.1_2.5m_bio_17.tif","../worldclim/cur_bio/wc2.1_2.5m_bio_18.tif",
                      "../worldclim/cur_bio/wc2.1_2.5m_bio_19.tif"))

f_245 <- terra::rast("../worldclim/wc2.1_2.5m_bioc_MIROC6_ssp245_2021-2040.tif")
f_585 <- terra::rast("../worldclim/wc2.1_2.5m_bioc_MIROC6_ssp585_2021-2040.tif")


names(curr) <- nms
names(f_245) <- nms
names(f_585) <- nms


bioclims <- terra::crop(curr, jp) %>%
  terra::mask(jp)

f_245 <- terra::crop(f_245, jp) %>%
  terra::mask(jp)

f_585 <- terra::crop(f_585, jp) %>%
  terra::mask(jp)


v1 <- usdm::vifstep(bioclims, th = 7 )

bio_tr <- exclude(bioclims, v1)
f_245_tr <- exclude(f_245, v1)
f_585_tr <- exclude(f_585, v1)

mdl_rf <- SSDM::modelling(c("RF"), Occurrences = carp_2_thined, Env = raster::stack(bio_tr),
                           Xcol = "long", Ycol = "lat")

mdl_rf <- read_rds("jp_R/mdl_rf_final2.rds")
v1  <- read_rds("jp_R/v1_final2.rds")


mdl <- read_rds("/home/borshon/infinity/jp_geo/jp_R/mdl_rf_final2.rds")


eval <- as.data.frame(mdl@evaluation %>%
                        mutate_if(is.numeric, round, 3))
eval <- tibble(eval[1],eval[2],eval[3],eval[4],eval[5],
               eval[6],eval[7],eval[8])

knitr::kable(eval,format = "markdown")


bioc <- tibble(
  "Codes" = c("bio3","bio4","bio5","bio8","bio15","bio18","bio19"),
  "Descriptions" = c("Isothermality","Temperature seasonality","Maximum Temp of warmest month","Mean Temp of Wettest Quarter",
                     "Precipitation Seasonality","Precipitation of Warmest Quarter","Precipitation of Coldest Quarter")
)

knitr::kable(bioc, format="markdown")


mdl@variable.importance %>%
  tidyr::gather(Variable, Value) %>%
  mutate(Variable = forcats::fct_reorder(Variable, Value)) %>%
  ggplot(aes(x = Variable, y = Value)) +
  geom_col(fill = '#82B2C3') +
  coord_flip() +
  geom_text(aes(label = round(Value, 2), y = Value + 0.18),
            size = 3, vjust = 1, hjust = -0.5 ) +
  labs(title = 'Variable Importance',
       subtitle = 'Common Carp RF model') +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(size = 14, hjust = 0,
                                  margin = margin(b = 2), face = 'bold'),
        plot.subtitle = element_text(size = 12, hjust = 0,
                                     margin = margin(b = 5), face = 'italic'),
        axis.text = element_text(color = 'black'))



curr_p <- mdl@projection %>% rast()
f_2_p <- SSDM::project(mdl, raster::stack(f_245_tr))@projection %>% rast()
f_5_p <- SSDM::project(mdl, raster::stack(f_585_tr))@projection %>% rast()

curr_pr <- mask(curr_p, jp_w)
f_2_pr <- mask(f_2_p, jp_w)
f_5_pr <- mask(f_5_p, jp_w)


fc <- freq(curr_p)
f2 <- freq(f_2_p)
f5 <- freq(f_5_p)

area <- (15.07+19.55)/2

fc$area = fc$count*area
f2$area = f2$count*area
f5$area = f5$count*area

distance <- tibble("climate change scenarios" = c("current","2100 ssp245","2100 ssp585"),
                   "range (km2)" = c(fc$area[2],f2$area[2],f5$area[2]),
                   "range increase (km2)" = c(fc$area[2]-fc$area[2],
                                              f2$area[2]-fc$area[2],
                                              f5$area[2]-fc$area[2]),
                   "range increase (%)" = c(as.integer((fc$area[2]/fc$area[2]*100)-100),
                                            as.integer((f2$area[2]/fc$area[2]*100)-100),
                                            as.integer((f5$area[2]/fc$area[2]*100)-100)))

knitr::kable(distance, format="markdown")


plot(curr_p, type="interval", breaks=c(0,0.2,0.5,0.7,1),  col=c("lightblue1","#8fe2d5","#607a9f","#020356"), main="On Mainland",
     plg=list(title="Probability", title.cex = 0.8),mar=c(3.1, 3.1, 2.1, 9))

plot(curr_pr, type="interval", breaks=c(0,0.2,0.5,0.7,1), col=c("lightblue1","#8fe2d5","#607a9f","#020356"), main="In Rivers and Streams",plg=list(title="Probability", title.cex = 0.8),mar=c(3.1, 3.1, 2.1, 9))

plot(f_2_p, type="interval", breaks=c(0,0.2,0.5,0.7,1),  col=c("lightblue1","#8fe2d5","#607a9f","#020356"), main="On Mainland",plg=list(title="Probability", title.cex = 0.8),mar=c(3.1, 3.1, 2.1, 9))

plot(f_2_pr, type="interval", breaks=c(0,0.2,0.5,0.7,1), col=c("lightblue1","#8fe2d5","#607a9f","#020356"), main="In Rivers and Streams",plg=list(title="Probability", title.cex = 0.8),mar=c(3.1, 3.1, 2.1, 9))

plot(f_5_p, type="interval", breaks=c(0,0.2,0.5,0.7,1),  col=c("lightblue1","#8fe2d5","#607a9f","#020356"), main="On Mainland",plg=list(title="Probability", title.cex = 0.8),mar=c(3.1, 3.1, 2.1, 9))

plot(f_5_pr, type="interval", breaks=c(0,0.2,0.5,0.7,1), col=c("lightblue1","#8fe2d5","#607a9f","#020356"), main="In Rivers and Streams",plg=list(title="Probability", title.cex = 0.8),mar=c(3.1, 3.1, 2.1, 9))










































