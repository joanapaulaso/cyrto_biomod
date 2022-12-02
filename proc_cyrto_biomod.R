# 01 - Call libraries ----
library(raster)
library(ENMTML)
library(biomod2)
library(rgdal)
library(sdmpredictors)
library(maptools)
library(usdm)
library(ecospat)
library(CoordinateCleaner)
library(rgbif)
library(spocc)
library(spThin)
library(dplyr)
library(dismo)
library(gridExtra)
library(tidyr)
library(corrplot)
library(geometry)
library(sp)
library(ggfortify)
library(ggplot2)

# 02 - Set wd ----
setwd("D:/SDM/cyrto22/known_distrib")
getwd()

# 03 - Preprocess ----
# Load environmental variables
# WorldClim
bio1 <- raster("./HIGHRES/wc2.1_30s_bio_1.tif")
bio2 <- raster("./HIGHRES/wc2.1_30s_bio_2.tif")
bio3 <- raster("./HIGHRES/wc2.1_30s_bio_3.tif")
bio4 <- raster("./HIGHRES/wc2.1_30s_bio_4.tif")
bio5 <- raster("./HIGHRES/wc2.1_30s_bio_5.tif")
bio6 <- raster("./HIGHRES/wc2.1_30s_bio_6.tif")
bio7 <- raster("./HIGHRES/wc2.1_30s_bio_7.tif")
bio8 <- raster("./HIGHRES/wc2.1_30s_bio_8.tif")
bio9 <- raster("./HIGHRES/wc2.1_30s_bio_9.tif")
bio10 <- raster("./HIGHRES/wc2.1_30s_bio_10.tif")
bio11 <- raster("./HIGHRES/wc2.1_30s_bio_11.tif")
bio12 <- raster("./HIGHRES/wc2.1_30s_bio_12.tif")
bio13 <- raster("./HIGHRES/wc2.1_30s_bio_13.tif")
bio14 <- raster("./HIGHRES/wc2.1_30s_bio_14.tif")
bio15 <- raster("./HIGHRES/wc2.1_30s_bio_15.tif")
bio16 <- raster("./HIGHRES/wc2.1_30s_bio_16.tif")
bio17 <- raster("./HIGHRES/wc2.1_30s_bio_17.tif")
bio18 <- raster("./HIGHRES/wc2.1_30s_bio_18.tif")
bio19 <- raster("./HIGHRES/wc2.1_30s_bio_19.tif")
prec <- raster("./HIGHRES/PREC_MEAN.tif")
srad <- raster("./HIGHRES/SRAD_MEAN.tif")
tavg <- raster("./HIGHRES/TAVG_MEAN.tif")
vapr <- raster("./HIGHRES/VAPR_MEAN.tif")
wind <- raster("./HIGHRES/WIND_MEAN.tif")


# Cloud cover MODCF
cc1 <- raster("./HIGHRES/MODCF_intraannualSD.tif")
cc2 <- raster("./HIGHRES/MODCF_monthlymean_01.tif")
cc3 <- raster("./HIGHRES/MODCF_monthlymean_07.tif")

# NDVI 2003-2007
ndvi1 <- raster("./HIGHRES/NDVI_CUM_2003-2007.tif")
ndvi2 <- raster("./HIGHRES/NDVI_MIN_2003-2007.tif")
ndvi3 <- raster("./HIGHRES/NDVI_SEAS_2003-2007.tif")

# SoilGrids
soil1 <- raster("./HIGHRES/CEC.tif") # Cation Exchange Capacity of the soil
soil2 <- raster("./HIGHRES/CFVO.tif") # Volumetric fraction of coarse fragments (> 2 mm)
soil3 <- raster("./HIGHRES/CLAY.tif")
soil4 <- raster("./HIGHRES/NITROGEN.tif")
soil5 <- raster("./HIGHRES/OCD.tif") # Organic carbon density
soil6 <- raster("./HIGHRES/phh2o.tif")
soil7 <- raster("./HIGHRES/SAND.tif")
soil8 <- raster("./HIGHRES/SILT.tif")
soil9 <- raster("./HIGHRES/SOC.tif") # Soil organic carbon content in the fine earth fraction

# Topography
topo1 <- raster("./HIGHRES/AI.tif")
topo2 <- raster("./HIGHRES/aspectcosine.tif")
topo3 <- raster("./HIGHRES/aspectsine.tif")
topo4 <- raster("./HIGHRES/CTI.tif")
topo5 <- raster("./HIGHRES/elevation.tif")
topo6 <- raster("./HIGHRES/HLI.tif")
topo7 <- raster("./HIGHRES/slope.tif")
topo8 <- raster("./HIGHRES/tpi.tif")
topo9 <- raster("./HIGHRES/vrm.tif")

# Create variables array
biovars = c(
  bio1, 
  bio2, 
  bio3, 
  bio4, 
  bio5, 
  bio6, 
  bio7, 
  bio8, 
  bio9, 
  bio10, 
  bio11, 
  bio12, 
  bio13, 
  bio14, 
  bio15, 
  bio16, 
  bio17, 
  bio18, 
  bio19,
  prec,
  srad,
  tavg,
  vapr,
  wind,
  cc1, 
  cc2, 
  cc3, 
  ndvi1, 
  ndvi2, 
  ndvi3, 
  soil1, 
  soil2, 
  soil3, 
  soil4, 
  soil5, 
  soil6, 
  soil7, 
  soil8, 
  soil9, 
  topo1, 
  topo2, 
  topo3, 
  topo4, 
  topo5, 
  topo6, 
  topo7, 
  topo8, 
  topo9
)

# Adjust resolutions
new_res <- bio16
bio1_R <- resample(bio1, new_res, method = 'bilinear')
bio2_R <- resample(bio2, new_res, method = 'bilinear')
bio3_R <- resample(bio3, new_res, method = 'bilinear')
bio4_R <- resample(bio4, new_res, method = 'bilinear')
bio5_R <- resample(bio5, new_res, method = 'bilinear')
bio6_R <- resample(bio6, new_res, method = 'bilinear')
bio7_R <- resample(bio7, new_res, method = 'bilinear')
bio8_R <- resample(bio8, new_res, method = 'bilinear')
bio9_R <- resample(bio9, new_res, method = 'bilinear')
bio10_R <- resample(bio10, new_res, method = 'bilinear')
bio11_R <- resample(bio11, new_res, method = 'bilinear')
bio12_R <- resample(bio12, new_res, method = 'bilinear')
bio13_R <- resample(bio13, new_res, method = 'bilinear')
bio14_R <- resample(bio14, new_res, method = 'bilinear')
bio15_R <- resample(bio15, new_res, method = 'bilinear')
bio16_R <- resample(bio16, new_res, method = 'bilinear')
bio17_R <- resample(bio17, new_res, method = 'bilinear')
bio18_R <- resample(bio18, new_res, method = 'bilinear')
bio19_R <- resample(bio19, new_res, method = 'bilinear')
prec_R <- resample(prec, new_res, method = 'bilinear')
srad_R <- resample(srad, new_res, method = 'bilinear')
tavg_R <- resample(tavg, new_res, method = 'bilinear')
vapr_R <- resample(vapr, new_res, method = 'bilinear')
wind_R <- resample(wind, new_res, method = 'bilinear')
cc1_R <- resample(cc1, new_res, method = 'bilinear')
cc2_R <- resample(cc2, new_res, method = 'bilinear')
cc3_R <- resample(cc3, new_res, method = 'bilinear')
ndvi1_R <- resample(ndvi1, new_res, method = 'bilinear')
ndvi2_R <- resample(ndvi2, new_res, method = 'bilinear')
ndvi3_R <- resample(ndvi3, new_res, method = 'bilinear')
soil1_R <- resample(soil1, new_res, method = 'bilinear')
soil2_R <- resample(soil2, new_res, method = 'bilinear')
soil3_R <- resample(soil3, new_res, method = 'bilinear')
soil4_R <- resample(soil4, new_res, method = 'bilinear')
soil5_R <- resample(soil5, new_res, method = 'bilinear')
soil6_R <- resample(soil6, new_res, method = 'bilinear')
soil7_R <- resample(soil7, new_res, method = 'bilinear')
soil8_R <- resample(soil8, new_res, method = 'bilinear')
soil9_R <- resample(soil9, new_res, method = 'bilinear')
topo1_R <- resample(topo1, new_res, method = 'bilinear')
topo2_R <- resample(topo2, new_res, method = 'bilinear')
topo3_R <- resample(topo3, new_res, method = 'bilinear')
topo4_R <- resample(topo4, new_res, method = 'bilinear')
topo5_R <- resample(topo5, new_res, method = 'bilinear')
topo6_R <- resample(topo6, new_res, method = 'bilinear')
topo7_R <- resample(topo7, new_res, method = 'bilinear')
topo8_R <- resample(topo8, new_res, method = 'bilinear')
topo9_R <- resample(topo9, new_res, method = 'bilinear')

# Crop variables (Brazil extent)
Brasil <- getData('GADM', country= 'Brazil', level=1, path = "./VAR")
bio1_R_BR <- mask(crop(bio1_R, Brasil), Brasil)
bio2_R_BR <- mask(crop(bio2_R, Brasil), Brasil)
bio3_R_BR <- mask(crop(bio3_R, Brasil), Brasil)
bio4_R_BR <- mask(crop(bio4_R, Brasil), Brasil)
bio5_R_BR <- mask(crop(bio5_R, Brasil), Brasil)
bio6_R_BR <- mask(crop(bio6_R, Brasil), Brasil)
bio7_R_BR <- mask(crop(bio7_R, Brasil), Brasil)
bio8_R_BR <- mask(crop(bio8_R, Brasil), Brasil)
bio9_R_BR <- mask(crop(bio9_R, Brasil), Brasil)
bio10_R_BR <- mask(crop(bio10_R, Brasil), Brasil)
bio11_R_BR <- mask(crop(bio11_R, Brasil), Brasil)
bio12_R_BR <- mask(crop(bio12_R, Brasil), Brasil)
bio13_R_BR <- mask(crop(bio13_R, Brasil), Brasil)
bio14_R_BR <- mask(crop(bio14_R, Brasil), Brasil)
bio15_R_BR <- mask(crop(bio15_R, Brasil), Brasil)
bio16_R_BR <- mask(crop(bio16_R, Brasil), Brasil)
bio17_R_BR <- mask(crop(bio17_R, Brasil), Brasil)
bio18_R_BR <- mask(crop(bio18_R, Brasil), Brasil)
bio19_R_BR <- mask(crop(bio19_R, Brasil), Brasil)
prec_R_BR <- mask(crop(prec_R, Brasil), Brasil)
srad_R_BR <- mask(crop(srad_R, Brasil), Brasil)
tavg_R_BR <- mask(crop(tavg_R, Brasil), Brasil)
vapr_R_BR <- mask(crop(vapr_R, Brasil), Brasil)
wind_R_BR <- mask(crop(wind_R, Brasil), Brasil)
cc1_R_BR <- mask(crop(cc1_R, Brasil), Brasil)
cc2_R_BR <- mask(crop(cc2_R, Brasil), Brasil)
cc3_R_BR <- mask(crop(cc3_R, Brasil), Brasil)
ndvi1_R_BR <- mask(crop(ndvi1_R, Brasil), Brasil)
ndvi2_R_BR <- mask(crop(ndvi2_R, Brasil), Brasil)
ndvi3_R_BR <- mask(crop(ndvi3_R, Brasil), Brasil)
soil1_R_BR <- mask(crop(soil1_R, Brasil), Brasil)
soil2_R_BR <- mask(crop(soil2_R, Brasil), Brasil)
soil3_R_BR <- mask(crop(soil3_R, Brasil), Brasil)
soil4_R_BR <- mask(crop(soil4_R, Brasil), Brasil)
soil5_R_BR <- mask(crop(soil5_R, Brasil), Brasil)
soil6_R_BR <- mask(crop(soil6_R, Brasil), Brasil)
soil7_R_BR <- mask(crop(soil7_R, Brasil), Brasil)
soil8_R_BR <- mask(crop(soil8_R, Brasil), Brasil)
soil9_R_BR <- mask(crop(soil9_R, Brasil), Brasil)
topo1_R_BR <- mask(crop(topo1_R, Brasil), Brasil)
topo2_R_BR <- mask(crop(topo2_R, Brasil), Brasil)
topo3_R_BR <- mask(crop(topo3_R, Brasil), Brasil)
topo4_R_BR <- mask(crop(topo4_R, Brasil), Brasil)
topo5_R_BR <- mask(crop(topo5_R, Brasil), Brasil)
topo6_R_BR <- mask(crop(topo6_R, Brasil), Brasil)
topo7_R_BR <- mask(crop(topo7_R, Brasil), Brasil)
topo8_R_BR <- mask(crop(topo8_R, Brasil), Brasil)
topo9_R_BR <- mask(crop(topo9_R, Brasil), Brasil)

# Save files with cropped raster
rasterDir = "./VAR/BrasilMask"
writeRaster(bio1_R_BR, filename = file.path(rasterDir, "bio1_R_BR.tif"), format = "GTiff")
writeRaster(bio2_R_BR, filename = file.path(rasterDir, "bio2_R_BR.tif"), format = "GTiff")
writeRaster(bio3_R_BR, filename = file.path(rasterDir, "bio3_R_BR.tif"), format = "GTiff")
writeRaster(bio4_R_BR, filename = file.path(rasterDir, "bio4_R_BR.tif"), format = "GTiff")
writeRaster(bio5_R_BR, filename = file.path(rasterDir, "bio5_R_BR.tif"), format = "GTiff")
writeRaster(bio6_R_BR, filename = file.path(rasterDir, "bio6_R_BR.tif"), format = "GTiff")
writeRaster(bio7_R_BR, filename = file.path(rasterDir, "bio7_R_BR.tif"), format = "GTiff")
writeRaster(bio8_R_BR, filename = file.path(rasterDir, "bio8_R_BR.tif"), format = "GTiff")
writeRaster(bio9_R_BR, filename = file.path(rasterDir, "bio9_R_BR.tif"), format = "GTiff")
writeRaster(bio10_R_BR, filename = file.path(rasterDir, "bio10_R_BR.tif"), format = "GTiff")
writeRaster(bio11_R_BR, filename = file.path(rasterDir, "bio11_R_BR.tif"), format = "GTiff")
writeRaster(bio12_R_BR, filename = file.path(rasterDir, "bio12_R_BR.tif"), format = "GTiff")
writeRaster(bio13_R_BR, filename = file.path(rasterDir, "bio13_R_BR.tif"), format = "GTiff")
writeRaster(bio14_R_BR, filename = file.path(rasterDir, "bio14_R_BR.tif"), format = "GTiff")
writeRaster(bio15_R_BR, filename = file.path(rasterDir, "bio15_R_BR.tif"), format = "GTiff")
writeRaster(bio16_R_BR, filename = file.path(rasterDir, "bio16_R_BR.tif"), format = "GTiff")
writeRaster(bio17_R_BR, filename = file.path(rasterDir, "bio17_R_BR.tif"), format = "GTiff")
writeRaster(bio18_R_BR, filename = file.path(rasterDir, "bio18_R_BR.tif"), format = "GTiff")
writeRaster(bio19_R_BR, filename = file.path(rasterDir, "bio19_R_BR.tif"), format = "GTiff")
writeRaster(prec_R_BR, filename = file.path(rasterDir, "prec_R_BR.tif"), format = "GTiff")
writeRaster(srad_R_BR, filename = file.path(rasterDir, "srad_R_BR.tif"), format = "GTiff")
writeRaster(tavg_R_BR, filename = file.path(rasterDir, "tavg_R_BR.tif"), format = "GTiff")
writeRaster(vapr_R_BR, filename = file.path(rasterDir, "vapr_R_BR.tif"), format = "GTiff")
writeRaster(wind_R_BR, filename = file.path(rasterDir, "wind_R_BR.tif"), format = "GTiff")
writeRaster(cc1_R_BR, filename = file.path(rasterDir, "cc1_R_BR.tif"), format = "GTiff")
writeRaster(cc2_R_BR, filename = file.path(rasterDir, "cc2_R_BR.tif"), format = "GTiff")
writeRaster(cc3_R_BR, filename = file.path(rasterDir, "cc3_R_BR.tif"), format = "GTiff")
writeRaster(ndvi1_R_BR, filename = file.path(rasterDir, "ndvi1_R_BR.tif"), format = "GTiff")
writeRaster(ndvi2_R_BR, filename = file.path(rasterDir, "ndvi2_R_BR.tif"), format = "GTiff")
writeRaster(ndvi3_R_BR, filename = file.path(rasterDir, "ndvi3_R_BR.tif"), format = "GTiff")
writeRaster(soil1_R_BR, filename = file.path(rasterDir, "soil1_R_BR.tif"), format = "GTiff")
writeRaster(soil2_R_BR, filename = file.path(rasterDir, "soil2_R_BR.tif"), format = "GTiff")
writeRaster(soil3_R_BR, filename = file.path(rasterDir, "soil3_R_BR.tif"), format = "GTiff")
writeRaster(soil4_R_BR, filename = file.path(rasterDir, "soil4_R_BR.tif"), format = "GTiff")
writeRaster(soil5_R_BR, filename = file.path(rasterDir, "soil5_R_BR.tif"), format = "GTiff")
writeRaster(soil6_R_BR, filename = file.path(rasterDir, "soil6_R_BR.tif"), format = "GTiff")
writeRaster(soil7_R_BR, filename = file.path(rasterDir, "soil7_R_BR.tif"), format = "GTiff")
writeRaster(soil8_R_BR, filename = file.path(rasterDir, "soil8_R_BR.tif"), format = "GTiff")
writeRaster(soil9_R_BR, filename = file.path(rasterDir, "soil9_R_BR.tif"), format = "GTiff")
writeRaster(topo1_R_BR, filename = file.path(rasterDir, "topo1_R_BR.tif"), format = "GTiff")
writeRaster(topo2_R_BR, filename = file.path(rasterDir, "topo2_R_BR.tif"), format = "GTiff")
writeRaster(topo3_R_BR, filename = file.path(rasterDir, "topo3_R_BR.tif"), format = "GTiff")
writeRaster(topo4_R_BR, filename = file.path(rasterDir, "topo4_R_BR.tif"), format = "GTiff")
writeRaster(topo5_R_BR, filename = file.path(rasterDir, "topo5_R_BR.tif"), format = "GTiff")
writeRaster(topo6_R_BR, filename = file.path(rasterDir, "topo6_R_BR.tif"), format = "GTiff")
writeRaster(topo7_R_BR, filename = file.path(rasterDir, "topo7_R_BR.tif"), format = "GTiff")
writeRaster(topo8_R_BR, filename = file.path(rasterDir, "topo8_R_BR.tif"), format = "GTiff")
writeRaster(topo9_R_BR, filename = file.path(rasterDir, "topo9_R_BR.tif"), format = "GTiff")

# 04 - Reload cropped files ----
bio1_R_BR <- raster("./VAR/BrasilMask/bio1_R_BR.tif")
bio2_R_BR <- raster("./VAR/BrasilMask/bio2_R_BR.tif")
bio3_R_BR <- raster("./VAR/BrasilMask/bio3_R_BR.tif")
bio4_R_BR <- raster("./VAR/BrasilMask/bio4_R_BR.tif")
bio5_R_BR <- raster("./VAR/BrasilMask/bio5_R_BR.tif")
bio6_R_BR <- raster("./VAR/BrasilMask/bio6_R_BR.tif")
bio7_R_BR <- raster("./VAR/BrasilMask/bio7_R_BR.tif")
bio8_R_BR <- raster("./VAR/BrasilMask/bio8_R_BR.tif")
bio9_R_BR <- raster("./VAR/BrasilMask/bio9_R_BR.tif")
bio10_R_BR <- raster("./VAR/BrasilMask/bio10_R_BR.tif")
bio11_R_BR <- raster("./VAR/BrasilMask/bio11_R_BR.tif")
bio12_R_BR <- raster("./VAR/BrasilMask/bio12_R_BR.tif")
bio13_R_BR <- raster("./VAR/BrasilMask/bio13_R_BR.tif")
bio14_R_BR <- raster("./VAR/BrasilMask/bio14_R_BR.tif")
bio15_R_BR <- raster("./VAR/BrasilMask/bio15_R_BR.tif")
bio16_R_BR <- raster("./VAR/BrasilMask/bio16_R_BR.tif")
bio17_R_BR <- raster("./VAR/BrasilMask/bio17_R_BR.tif")
bio18_R_BR <- raster("./VAR/BrasilMask/bio18_R_BR.tif")
bio19_R_BR <- raster("./VAR/BrasilMask/bio19_R_BR.tif")
prec_R_BR <- raster("./VAR/BrasilMask/prec_R_BR.tif")
srad_R_BR <- raster("./VAR/BrasilMask/srad_R_BR.tif")
tavg_R_BR <- raster("./VAR/BrasilMask/tavg_R_BR.tif")
vapr_R_BR <- raster("./VAR/BrasilMask/vapr_R_BR.tif")
wind_R_BR <- raster("./VAR/BrasilMask/wind_R_BR.tif")
prec_R_BR <- raster("./VAR/BrasilMask/prec_R_BR.tif")
srad_R_BR <- raster("./VAR/BrasilMask/srad_R_BR.tif")
tavg_R_BR <- raster("./VAR/BrasilMask/tavg_R_BR.tif")
vapr_R_BR <- raster("./VAR/BrasilMask/vapr_R_BR.tif")
wind_R_BR <- raster("./VAR/BrasilMask/wind_R_BR.tif")
cc1_R_BR <- raster("./VAR/BrasilMask/cc1_R_BR.tif")
cc2_R_BR <- raster("./VAR/BrasilMask/cc2_R_BR.tif")
cc3_R_BR <- raster("./VAR/BrasilMask/cc3_R_BR.tif")
ndvi1_R_BR <- raster("./VAR/BrasilMask/ndvi1_R_BR.tif")
ndvi2_R_BR <- raster("./VAR/BrasilMask/ndvi2_R_BR.tif")
ndvi3_R_BR <- raster("./VAR/BrasilMask/ndvi3_R_BR.tif")
soil1_R_BR <- raster("./VAR/BrasilMask/soil1_R_BR.tif")
soil2_R_BR <- raster("./VAR/BrasilMask/soil2_R_BR.tif")
soil3_R_BR <- raster("./VAR/BrasilMask/soil3_R_BR.tif")
soil4_R_BR <- raster("./VAR/BrasilMask/soil4_R_BR.tif")
soil5_R_BR <- raster("./VAR/BrasilMask/soil5_R_BR.tif")
soil6_R_BR <- raster("./VAR/BrasilMask/soil6_R_BR.tif")
soil7_R_BR <- raster("./VAR/BrasilMask/soil7_R_BR.tif")
soil8_R_BR <- raster("./VAR/BrasilMask/soil8_R_BR.tif")
soil9_R_BR <- raster("./VAR/BrasilMask/soil9_R_BR.tif")
topo1_R_BR <- raster("./VAR/BrasilMask/topo1_R_BR.tif")
topo2_R_BR <- raster("./VAR/BrasilMask/topo2_R_BR.tif")
topo3_R_BR <- raster("./VAR/BrasilMask/topo3_R_BR.tif")
topo4_R_BR <- raster("./VAR/BrasilMask/topo4_R_BR.tif")
topo5_R_BR <- raster("./VAR/BrasilMask/topo5_R_BR.tif")
topo6_R_BR <- raster("./VAR/BrasilMask/topo6_R_BR.tif")
topo7_R_BR <- raster("./VAR/BrasilMask/topo7_R_BR.tif")
topo8_R_BR <- raster("./VAR/BrasilMask/topo8_R_BR.tif")
topo9_R_BR <- raster("./VAR/BrasilMask/topo9_R_BR.tif")

# 05 - Create biostacks ----
biostack_WC <- stack(
  bio1_R_BR,
  bio2_R_BR,
  bio3_R_BR,
  bio4_R_BR,
  bio5_R_BR,
  bio6_R_BR,
  bio7_R_BR,
  bio8_R_BR,
  bio9_R_BR,
  bio10_R_BR,
  bio11_R_BR,
  bio12_R_BR,
  bio13_R_BR,
  bio14_R_BR,
  bio15_R_BR,
  bio16_R_BR,
  bio17_R_BR,
  bio18_R_BR,
  bio19_R_BR,
  prec_R_BR,
  srad_R_BR,
  tavg_R_BR,
  vapr_R_BR,
  wind_R_BR
)

biostack_Topographic <- stack(
  topo1_R_BR,
  topo2_R_BR,
  topo3_R_BR,
  topo4_R_BR,
  topo5_R_BR,
  topo6_R_BR,
  topo7_R_BR,
  topo8_R_BR,
  topo9_R_BR
)

biostack_Soil <- stack(
  soil1_R_BR,
  soil2_R_BR,
  soil3_R_BR,
  soil4_R_BR,
  soil5_R_BR,
  soil6_R_BR,
  soil7_R_BR,
  soil8_R_BR,
  soil9_R_BR
)

biostack_CloudCover <- stack(
  cc1_R_BR,
  cc2_R_BR,
  cc3_R_BR
)

biostack_NDVI <- stack(
  ndvi1_R_BR,
  ndvi2_R_BR,
  ndvi3_R_BR
)

biostack_allvar <- stack(
  bio1_R_BR,
  bio2_R_BR,
  bio3_R_BR,
  bio4_R_BR,
  bio5_R_BR,
  bio6_R_BR,
  bio7_R_BR,
  bio8_R_BR,
  bio9_R_BR,
  bio10_R_BR,
  bio11_R_BR,
  bio12_R_BR,
  bio13_R_BR,
  bio14_R_BR,
  bio15_R_BR,
  bio16_R_BR,
  bio17_R_BR,
  bio18_R_BR,
  bio19_R_BR,
  prec_R_BR,
  srad_R_BR,
  tavg_R_BR,
  vapr_R_BR,
  wind_R_BR,
  topo1_R_BR,
  topo2_R_BR,
  topo3_R_BR,
  topo4_R_BR,
  topo5_R_BR,
  topo6_R_BR,
  topo7_R_BR,
  topo8_R_BR,
  topo9_R_BR,
  soil1_R_BR,
  soil2_R_BR,
  soil3_R_BR,
  soil4_R_BR,
  soil5_R_BR,
  soil6_R_BR,
  soil7_R_BR,
  soil8_R_BR,
  soil9_R_BR,
  cc1_R_BR,
  cc2_R_BR,
  cc3_R_BR,
  ndvi1_R_BR,
  ndvi2_R_BR,
  ndvi3_R_BR
)

# 06 - Get polygon ----
polygonDir <- "./POLYGON"
cyrto_polygon <- readOGR(polygonDir, "clipped_cg")
plot(prec_R_BR, main="Bio 1 Brasil")
plot(cyrto_polygon, add = T)

# 07 - Crop variables (cyrto_polygon extent) ----
bio1_cyrto <- mask(crop(bio1_R_BR, cyrto_polygon), cyrto_polygon)
bio2_cyrto <- mask(crop(bio2_R_BR, cyrto_polygon), cyrto_polygon)
bio3_cyrto <- mask(crop(bio3_R_BR, cyrto_polygon), cyrto_polygon)
bio4_cyrto <- mask(crop(bio4_R_BR, cyrto_polygon), cyrto_polygon)
bio5_cyrto <- mask(crop(bio5_R_BR, cyrto_polygon), cyrto_polygon)
bio6_cyrto <- mask(crop(bio6_R_BR, cyrto_polygon), cyrto_polygon)
bio7_cyrto <- mask(crop(bio7_R_BR, cyrto_polygon), cyrto_polygon)
bio8_cyrto <- mask(crop(bio8_R_BR, cyrto_polygon), cyrto_polygon)
bio9_cyrto <- mask(crop(bio9_R_BR, cyrto_polygon), cyrto_polygon)
bio10_cyrto <- mask(crop(bio10_R_BR, cyrto_polygon), cyrto_polygon)
bio11_cyrto <- mask(crop(bio11_R_BR, cyrto_polygon), cyrto_polygon)
bio12_cyrto <- mask(crop(bio12_R_BR, cyrto_polygon), cyrto_polygon)
bio13_cyrto <- mask(crop(bio13_R_BR, cyrto_polygon), cyrto_polygon)
bio14_cyrto <- mask(crop(bio14_R_BR, cyrto_polygon), cyrto_polygon)
bio15_cyrto <- mask(crop(bio15_R_BR, cyrto_polygon), cyrto_polygon)
bio16_cyrto <- mask(crop(bio16_R_BR, cyrto_polygon), cyrto_polygon)
bio17_cyrto <- mask(crop(bio17_R_BR, cyrto_polygon), cyrto_polygon)
bio18_cyrto <- mask(crop(bio18_R_BR, cyrto_polygon), cyrto_polygon)
bio19_cyrto <- mask(crop(bio19_R_BR, cyrto_polygon), cyrto_polygon)
prec_cyrto <- mask(crop(prec_R_BR, cyrto_polygon), cyrto_polygon)
srad_cyrto <- mask(crop(srad_R_BR, cyrto_polygon), cyrto_polygon)
tavg_cyrto <- mask(crop(tavg_R_BR, cyrto_polygon), cyrto_polygon)
vapr_cyrto <- mask(crop(vapr_R_BR, cyrto_polygon), cyrto_polygon)
wind_cyrto <- mask(crop(wind_R_BR, cyrto_polygon), cyrto_polygon)
cc1_cyrto <- mask(crop(cc1_R_BR, cyrto_polygon), cyrto_polygon)
cc2_cyrto <- mask(crop(cc2_R_BR, cyrto_polygon), cyrto_polygon)
cc3_cyrto <- mask(crop(cc3_R_BR, cyrto_polygon), cyrto_polygon)
ndvi1_cyrto <- mask(crop(ndvi1_R_BR, cyrto_polygon), cyrto_polygon)
ndvi2_cyrto <- mask(crop(ndvi2_R_BR, cyrto_polygon), cyrto_polygon)
ndvi3_cyrto <- mask(crop(ndvi3_R_BR, cyrto_polygon), cyrto_polygon)
soil1_cyrto <- mask(crop(soil1_R_BR, cyrto_polygon), cyrto_polygon)
soil2_cyrto <- mask(crop(soil2_R_BR, cyrto_polygon), cyrto_polygon)
soil3_cyrto <- mask(crop(soil3_R_BR, cyrto_polygon), cyrto_polygon)
soil4_cyrto <- mask(crop(soil4_R_BR, cyrto_polygon), cyrto_polygon)
soil5_cyrto <- mask(crop(soil5_R_BR, cyrto_polygon), cyrto_polygon)
soil6_cyrto <- mask(crop(soil6_R_BR, cyrto_polygon), cyrto_polygon)
soil7_cyrto <- mask(crop(soil7_R_BR, cyrto_polygon), cyrto_polygon)
soil8_cyrto <- mask(crop(soil8_R_BR, cyrto_polygon), cyrto_polygon)
soil9_cyrto <- mask(crop(soil9_R_BR, cyrto_polygon), cyrto_polygon)
topo1_cyrto <- mask(crop(topo1_R_BR, cyrto_polygon), cyrto_polygon)
topo2_cyrto <- mask(crop(topo2_R_BR, cyrto_polygon), cyrto_polygon)
topo3_cyrto <- mask(crop(topo3_R_BR, cyrto_polygon), cyrto_polygon)
topo4_cyrto <- mask(crop(topo4_R_BR, cyrto_polygon), cyrto_polygon)
topo5_cyrto <- mask(crop(topo5_R_BR, cyrto_polygon), cyrto_polygon)
topo6_cyrto <- mask(crop(topo6_R_BR, cyrto_polygon), cyrto_polygon)
topo7_cyrto <- mask(crop(topo7_R_BR, cyrto_polygon), cyrto_polygon)
topo8_cyrto <- mask(crop(topo8_R_BR, cyrto_polygon), cyrto_polygon)
topo9_cyrto <- mask(crop(topo9_R_BR, cyrto_polygon), cyrto_polygon)

biostack_allvar_buffer <- stack(
  bio1_cyrto,
  bio2_cyrto,
  bio3_cyrto,
  bio4_cyrto,
  bio5_cyrto,
  bio6_cyrto,
  bio7_cyrto,
  bio8_cyrto,
  bio9_cyrto,
  bio10_cyrto,
  bio11_cyrto,
  bio12_cyrto,
  bio13_cyrto,
  bio14_cyrto,
  bio15_cyrto,
  bio16_cyrto,
  bio17_cyrto,
  bio18_cyrto,
  bio19_cyrto,
  prec_cyrto,
  srad_cyrto,
  tavg_cyrto,
  vapr_cyrto,
  wind_cyrto,
  topo1_cyrto,
  topo2_cyrto,
  topo3_cyrto,
  topo4_cyrto,
  topo5_cyrto,
  topo6_cyrto,
  topo7_cyrto,
  topo8_cyrto,
  topo9_cyrto,
  soil1_cyrto,
  soil2_cyrto,
  soil3_cyrto,
  soil4_cyrto,
  soil5_cyrto,
  soil6_cyrto,
  soil7_cyrto,
  soil8_cyrto,
  soil9_cyrto,
  cc1_cyrto,
  cc2_cyrto,
  cc3_cyrto,
  ndvi1_cyrto,
  ndvi2_cyrto,
  ndvi3_cyrto
)

# 08 - Crop biostack (cyrto_polygon extent) ----
biostack_WC_buffer <- mask(crop(biostack_WC, cyrto_polygon), cyrto_polygon)
biostack_Topographic_buffer <- mask(crop(biostack_Topographic, cyrto_polygon), cyrto_polygon)
biostack_Soil_buffer <- mask(crop(biostack_Soil, cyrto_polygon), cyrto_polygon)
biostack_CloudCover_buffer <- mask(crop(biostack_CloudCover, cyrto_polygon), cyrto_polygon)
biostack_NDVI_buffer <- mask(crop(biostack_NDVI, cyrto_polygon), cyrto_polygon)

# 09 - Get and process distribution coordinates ----

cyrtopodium_matrix <- read.csv("cyrto_distrib.csv")
summary(cyrtopodium_matrix)
cyrtopodium_matrix <- cyrtopodium_matrix %>% tidyr::drop_na(decimalLongitude, decimalLatitude)
flags_spatial <- CoordinateCleaner::clean_coordinates(
  x = cyrtopodium_matrix, 
  species = "species",
  lon = "decimalLongitude", 
  lat = "decimalLatitude",
  tests = c("capitals",
            "centroids",
            "duplicates",
            "equal",
            "gbif",
            "institutions",
            "seas",
            "urban",
            "validity",
            "zeros"
  )
)
head(flags_spatial)
summary(flags_spatial)
cyrtopodium_matrix <- cyrtopodium_matrix %>% 
  dplyr::filter(flags_spatial$.summary == TRUE)

plot(cc3_cyrto, xlim=c(-20,30),ylim=c(10,60), main="Pontos de ocorrência X Bio13")
data(wrld_simpl)
plot(wrld_simpl, add = TRUE, border = 'darkgray', lwd = 1)
points(cyrtopodium_matrix$decimalLongitude, cyrtopodium_matrix$decimalLatitude, pch = 21, bg = 'blue', cex = 1)

dir.create("./thinned")
cyrtopodium_matrix_thin <- 
  thin( 
    loc.data = cyrtopodium_matrix, 
    lat.col = "decimalLatitude", long.col = "decimalLongitude", 
    spec.col = "species", 
    thin.par = 10,
    reps = 10, 
    locs.thinned.list.return = TRUE, 
    write.files = TRUE, 
    max.files = 2, 
    out.dir = "./thinned", out.base = "cyrtopodium", 
    write.log.file = TRUE,
    log.file = "./thinned/cyrtopodium_log.txt"
  )

plotThin(cyrtopodium_matrix_thin)
cyrtopodium_matrix_thin = read.csv("./thinned/cyrtopodium_thin1.csv", sep=",")
View(cyrtopodium_matrix_thin)
summary(cyrtopodium_matrix_thin)

plot(cc3_cyrto, xlim=c(-20,30),ylim=c(10,60), main="Pontos de ocorrência X CC3")
data(wrld_simpl)
plot(wrld_simpl, add = TRUE, border = 'darkgray', lwd = 1)
points(cyrtopodium_matrix$decimalLongitude, cyrtopodium_matrix$decimalLatitude, pch = 21, bg = 'blue', cex = 1)
points(cyrtopodium_matrix_thin$decimalLongitude, cyrtopodium_matrix_thin$decimalLatitude, pch = 21, bg = 'red', cex = 1)

# 10 - Correlations ----

# TESTAR A IMPORTÂNCIA DE VARIÁVEIS ANTES DA CORRELAÇÃO (OK)

set.seed(2152)

# Selected from VAR_IMP
biostack_varimpSelected_buffer <- stack(
  cc1_cyrto,
  cc3_cyrto,
  bio3_cyrto,
  ndvi3_cyrto,
  srad_cyrto,
  bio13_cyrto,
  ndvi1_cyrto,
  soil2_cyrto,
  bio4_cyrto,
  topo7_cyrto,
  bio17_cyrto,
  topo9_cyrto,
  bio14_cyrto,
  soil8_cyrto,
  bio19_cyrto,
  bio15_cyrto,
  soil4_cyrto,
  topo5_cyrto,
  ndvi2_cyrto,
  bio8_cyrto,
  soil9_cyrto,
  bio5_cyrto
)

backgr_varimpSelected <- randomPoints(biostack_varimpSelected_buffer, 10000)
absclim_varimpSelected <- data.frame(raster::extract(biostack_varimpSelected_buffer, backgr_varimpSelected))
absclim_varimpSelected_std <- data.frame(scale(absclim_varimpSelected))
M_varimpSelected <- cor(absclim_varimpSelected_std)
write.csv(M_varimpSelected, "cor_varimpSelected_24-12.csv")
corrplot(M_varimpSelected, upper = "square", lower = "shade", insig='blank', number.cex = 0.8)
vifcor(biostack_varimpSelected_buffer, th = 0.7) # vif_cor.txt
vifstep(biostack_varimpSelected_buffer, th = 10) # vif_cor.txt

# 11 - Stack variables selection from groups (plus NDVI and CloudCover) ----
biostack_varSelection_noCorr <- stack(
  cc1_cyrto,
  bio3_cyrto,
  ndvi3_cyrto,
  srad_cyrto,
  bio13_cyrto,
  soil2_cyrto,
  bio4_cyrto,
  topo7_cyrto,
  bio14_cyrto
)
plot(biostack_varSelection_noCorr)

# 12 - Correlations from varSelection_noCorr (with NDVI and CloudCover) ----
backgr_Selected <- randomPoints(biostack_varSelection_noCorr, 10000)
absclim_Selected <- data.frame(raster::extract(biostack_varSelection_noCorr, backgr_Selected))
absclim_Selected_std <- data.frame(scale(absclim_Selected))
M_Selected <- cor(absclim_Selected_std)
write.csv(M_Selected, "cor_Selected_for_PCA_24-12.csv")
corrplot(M_Selected, upper = "square", lower = "shade", insig='blank', number.cex = 0.8) # cor1.pdf
vifcor(biostack_varSelection_noCorr, th = 0.7) # vif_cor.txt
vifstep(biostack_varSelection_noCorr, th = 10) # vif_cor.txt

# Principal Component Analysis of selected data (as above)
write.csv(absclim_Selected_std, "selected_df_for_PCA_24-12.csv")
selected_df <- read.csv("selected_df_for_PCA_24-12.csv")
selected_df <- selected_df %>% tidyr::drop_na()
selected_df <- selected_df[-1]
selected_df_pca <- princomp(selected_df, cor = TRUE, scores = TRUE)
write.csv(selected_df_pca[['loadings']], "selected_df_pcaLoadings_for_PCA_24-12.csv")
write.csv(selected_df_pca[['scores']], "selected_df_pcaScores_for_PCA_24-12.csv")
write.csv(selected_df_pca[['sdev']], "selected_df_pcaSdev_for_PCA_24-12.csv")
plot(selected_df_pca)
biplot(selected_df_pca, pc.biplot = TRUE)

# Criteria: corr > vfi > pca

# 13 - Stack final variables selection ----
biostack_final <- stack(
  cc1_cyrto,
  bio3_cyrto,
  ndvi3_cyrto,
  srad_cyrto,
  soil2_cyrto,
  topo7_cyrto
)
plot(biostack_final)