#                             --- dsdim() PRACTISE ---

## load DSsim package
rm(list=ls())
library(DSsim)

## simulate study area
coords <- gaps <- list()
coords[[1]] <- list(data.frame(x = c(0,1000,1000,0,0), y = c(0,0,
                                                             1000,1000,0)))
gaps[[1]] <- list(data.frame(x = c(400,600,500,350,400), y = c(100,
                                                               250,600,120,100)))
region <- make.region(region.name = "study.area", units = "m",
                      coords = coords, gaps = gaps)
plot(region)

data(transects.shp)

shapefile.pathway <- "C:/Users/Matthew/Documents/SUMMER SCHOOL/MRDS_bootstrap_variance_estimation/dssim"
write.shapefile(transects.shp, paste(shapefile.pathway,"/transects_1",
                                     sep = ""))
parallel.design <- make.design(transect.type = "Line", 
                               design.details = c("Parallel","Systematic"), region = region, 
                               design.axis = 0, spacing = 100, plus.sampling =FALSE, 
                               path = shapefile.pathway)
pop.density <- make.density(region.obj = region, x.space = 10, 
                            y.space = 10, constant = 0.5) 
pop.density <- add.hotspot(pop.density, centre = c(50, 200), 
                           sigma = 100, amplitude = 0.1)
pop.density <- add.hotspot(pop.density, centre = c(500, 700), 
                           sigma = 900, amplitude = 0.05)
pop.density <- add.hotspot(pop.density, centre = c(300, 100), 
                           sigma = 100, amplitude = -0.15)

pop.description <- make.population.description(N = 100, 
                                               density.obj = pop.density, region = region, fixed.N = TRUE)

## detection function and survey results
scale.par = 25

for (i in 1:2) {
  set.seed(123456)
  detect <- make.detectability(key.function = "hn", scale.param = scale.par,
                             truncation = 30) 
  ddf.analyses <- make.ddf.analysis.list(dsmodel = list(~cds(key = "hn",
                                                             formula = ~1),~cds(key = "hr", formula = ~1)), method = "ds", 
                                         criteria = "AIC")
  my.simulation <- make.simulation(reps = 10, single.transect.set = TRUE,
                                   region.obj = region, design.obj = parallel.design, 
                                   population.description.obj = pop.description, 
                                   detectability.obj = detect, ddf.analyses.list = ddf.analyses)
  if (i == 1) {
    survey.results <- create.survey.results(my.simulation, dht.table = TRUE)
    data <- survey.results@ddf.data@ddf.dat
    obs.table <- survey.results@obs.table@obs.table
    objects.1 <- unique(data$object)
    scale.par <- scale.par - 10
  } else {
    survey.results.2 <- create.survey.results(my.simulation, dht.table = TRUE)
    data.2 <- survey.results.2@ddf.data@ddf.dat
    obs.table.2 <- survey.results.2@obs.table@obs.table
    data <- unique(rbind(data, data.2))
    obs.table <- unique(rbind(obs.table, obs.table.2))
    
    objects.2 <- unique(data.2$object)
    all.objects <- data$object
    det.1 <- as.numeric(all.objects %in% objects.1)
    det.2 <- as.numeric(all.objects %in% objects.2)
    
    data <- data[rep(seq_len(nrow(data)), each=2),]
    detected <- c(rbind(det.1,det.2))
    data <- cbind(data, detected)
    
    ## tidy up data
    data <- subset(data, select=-c(x,y))
    names(data)[names(data)=="transect.ID"] <- "Sample.Label"
    obs.table <- obs.table[order(obs.table$Sample.Label),]
    region.table <- survey.results@region.table@region.table
    sample.table <- survey.results@sample.table@sample.table
  }
}

## dht tables
tables <- list("data" = data,
               "region.table" = region.table,
               "sample.table" = sample.table,
               "obs.table" = obs.table)
#--

