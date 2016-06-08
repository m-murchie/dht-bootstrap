#               --- DOUBLE OBSERVER DATA SIM FUNCTION: VERSION 1 ---

double.obs.sim <- function(region.obj, design.obj, pop.description.obj, detect.obj.1,
                           detect.obj.2, ddf.analyses.list, seed = 123456, 
                           plot=FALSE) {
  set.seed(seed)
    
  ## simulation object and survery results for observer 1
  my.simulation <- make.simulation(reps = 10, single.transect.set = TRUE,
                                   double.observer = FALSE, region.obj, design.obj, pop.description.obj, 
                                   detect.obj.1, ddf.analyses.list)
  survey.results <- create.survey.results(my.simulation, dht.table = TRUE)
  data <- survey.results@ddf.data@ddf.dat
  obs.table <- survey.results@obs.table@obs.table
  objects.1 <- unique(data$object)
      
  ## as above, observer 2
  set.seed(seed)
    
  my.simulation.2 <- make.simulation(reps = 10, single.transect.set = TRUE,
                                     double.observer = FALSE, region.obj, design.obj, pop.description.obj, 
                                     detect.obj.2, ddf.analyses.list)
  survey.results.2 <- create.survey.results(my.simulation.2, dht.table = TRUE)
  data.2 <- survey.results.2@ddf.data@ddf.dat
  obs.table.2 <- survey.results.2@obs.table@obs.table
  data <- unique(rbind(data, data.2))
  obs.table <- unique(rbind(obs.table, obs.table.2))
    
  ## objects dectected by both observers
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
      
  ## dht tables
  tables <- list("data" = data,
                 "region.table" = region.table, 
                 "sample.table" = sample.table,
                 "obs.table" = obs.table)
  
  if (plot == TRUE) {
    plot(survey.results)
    plot(survey.results.2)
  }
      
  return(tables)
}