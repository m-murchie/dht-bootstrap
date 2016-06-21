#               --- DOUBLE OBSERVER DATA SIM FUNCTION: PRACTISE 2 ---

## load DSsim package and source double observer data sim function
rm(list=ls())
library(DSsim)
setwd("C:/Users/Matthew/Documents/SUMMER SCHOOL/MRDS_bootstrap_variance_estimation")
source("R/double.obs.sim.R")
source("R/bootstrap function.R")

## set directory
setwd("C:/Users/Matthew/Documents/SUMMER SCHOOL/MRDS_bootstrap_variance_estimation/DSsim Exercise")

## region object
region.shapefile <- read.shapefile("Region")
region <- make.region(region.name = "Survey Region", units = "m", 
                      shapefile = region.shapefile)

## density object
load("density.surface.robj")
pop.density <- make.density(region = region, density.surface = density.surface, 
                            x.space = 1000, y.space = 1000) 

## transects and survey design
design.path <-"C:/Users/Matthew/Documents/SUMMER SCHOOL/MRDS_bootstrap_variance_estimation/DSsim Exercise/Survey Transects/Zigzag Design"
subjective.design <- make.design(transect.type = "Line", 
                                 design.details = c("user specified"), 
                                 region = region, plus.sampling = FALSE, 
                                 path = design.path)

## population description 
pop.description <- make.population.description(region.obj = region, 
                                               density.obj = pop.density, 
                                               N = 1000, fixed.N = TRUE)

## detectability for each observer and analyses object
detect.1 <- make.detectability(key.function = "hn", scale.param = 3000, 
                               truncation = 1000)
detect.2 <- make.detectability(key.function = "hn", scale.param = 2800, 
                               truncation = 1000)
ddf.analyses <- make.ddf.analysis.list(
  dsmodel = list(~cds(key = "hn", formula = ~1),   #half-normal model
                 ~cds(key = "hr", formula = ~1)),  #hazard rate model
  method = "ds", criteria = "AIC")

#---

tables <- double.obs.sim(region, subjective.design, pop.description, detect.1, detect.2, 
                         ddf.analyses, plot=TRUE)


ddf.model <- ddf(method = 'io.fi', mrmodel=~glm(link='logit', formula=~distance), 
                 data = tables$data, meta.data=list(width=1000))
dht.results <- dht(ddf.model, tables$region.table, tables$sample.table, 
                   tables$obs.table)


boot.dht(tables, trunc = 1000, hist = TRUE)
