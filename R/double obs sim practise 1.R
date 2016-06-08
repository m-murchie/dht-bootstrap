#               --- DOUBLE OBSERVER DATA SIM FUNCTION: PRACTISE 1 ---

## load DSsim package and source double observer data sim function
rm(list=ls())
library(DSsim)
setwd("C:/Users/Matthew/Documents/SUMMER SCHOOL/MRDS_bootstrap_variance_estimation")
source("R/double.obs.sim.R")

## region object
coords <- gaps <- list()
coords[[1]] <- list(data.frame(x = c(0,1000,1000,0,0), y = c(0,0,
                                                             1000,1000,0)))
gaps[[1]] <- list(data.frame(x = c(400,600,500,350,400), y = c(100,
                                                               250,600,120,100)))
region <- make.region(region.name = "study.area", units = "m",
                      coords = coords, gaps = gaps)

## density object
pop.density <- make.density(region.obj = region, x.space = 10, 
                            y.space = 10, constant = 0.5) 
pop.density <- add.hotspot(pop.density, centre = c(50, 200), 
                           sigma = 100, amplitude = 0.1)
pop.density <- add.hotspot(pop.density, centre = c(500, 700), 
                           sigma = 900, amplitude = 0.05)
pop.density <- add.hotspot(pop.density, centre = c(300, 100), 
                           sigma = 100, amplitude = -0.15)

## transects and survey design
data(transects.shp)

shapefile.pathway <- "C:/Users/Matthew/Documents/SUMMER SCHOOL/MRDS_bootstrap_variance_estimation/dssim"
write.shapefile(transects.shp, paste(shapefile.pathway,"/transects_1",
                                     sep = ""))
parallel.design <- make.design(transect.type = "Line", 
                               design.details = c("Parallel","Systematic"), region = region, 
                               design.axis = 0, spacing = 100, plus.sampling =FALSE, 
                               path = shapefile.pathway)

## population description
pop.description <- make.population.description(N = 100, 
                                               density.obj = pop.density, region = region, fixed.N = TRUE)

## detectability for each observer and analyses object
detect.1 <- make.detectability(key.function = "hn", scale.param = 25,
                               truncation = 30)
detect.2 <- make.detectability(key.function = "hn", scale.param = 15,
                               truncation = 30)

ddf.analyses <- make.ddf.analysis.list(dsmodel = list(~cds(key = "hn",
                                                           formula = ~1),~cds(key = "hr", formula = ~1)), method = "ds", 
                                       criteria = "AIC")

#---

double.obs.sim(region, parallel.design, pop.description, detect.1, detect.2, 
               ddf.analyses, plot=TRUE)

