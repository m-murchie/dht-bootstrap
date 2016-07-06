#               --- DOUBLE OBSERVER DATA SIM FUNCTION: VERSION 1 ---

double.obs.sim.mod.2 <- function(region.obj, design.obj, pop.description.obj, 
                                 detect.obj.1, detect.obj.2, ddf.analyses.list, 
                                 seed = 123456, plot=FALSE) {
  
  set.seed(seed)
    
  ## simulation object and survey results for observer 1
  my.simulation <- make.simulation(reps = 1, single.transect.set = TRUE,
                                   double.observer = FALSE, region.obj, design.obj, 
                                   pop.description.obj, detect.obj.1, 
                                   ddf.analyses.list)
  survey.results <- create.survey.results(my.simulation, dht.table = TRUE)
  population <- survey.results@population@population
  data <- survey.results@ddf.data@ddf.dat
  objects.1 <- unique(data$object)
      
  ## first observer objects detected by second observer
  det.probs.seen <- mrds:::keyfct.hn(data$distance, detect.obj.2@scale.param)
  rand.numbs.seen <- runif(length(det.probs.seen))
  det.2.seen <- as.numeric(rand.numbs.seen < det.probs.seen)
  
  ## objects not detected by first observer, within truncation distance
  truncation <- survey.results@population@detectability@truncation
  poss.detect <- DSsim:::calc.poss.detect.dists(survey.results@population, 
                                                survey.results@transects, 
                                                truncation)
  obj.missed <- poss.detect[!(poss.detect$object %in% data$object),]
  obj.missed <- subset(obj.missed, select=-available)
  
  ## missed objects detected by second observer
  det.probs.miss <- DSsim:::hn.detect(obj.missed$distance, detect.obj.2)
  rand.numbs.miss <- runif(length(det.probs.miss))
  det.2.miss <- as.numeric(rand.numbs.miss < det.probs.miss)
  obj.missed <- obj.missed[as.logical(det.2.miss),]

  ## detections
  det.1 <- c(rep(1, length(data$object)), rep(0, length(obj.missed$object)))
  det.2 <- append(det.2.seen, rep(1, length(obj.missed$object)))
  data <- rbind(data, obj.missed)
  data <- data[rep(seq_len(nrow(data)), each=2),]
  detected <- c(rbind(det.1,det.2))
  
  ## tidy up data
  data <- subset(data, select=-c(x,y))
  observer <- rep(1:2, length(detected)/2)
  data <- cbind(data["object"], observer, detected, data[,2:3])
  names(data)[names(data)=="transect.ID"] <- "Sample.Label"
  
  
  ## add region column to data
  data["Region.Label"] <- rep("study.area", length(data$object))
  region.table <- survey.results@region.table@region.table
  sample.table <- survey.results@sample.table@sample.table
  
  ## dht tables
  tables <- list("true.Nc"      = length(poss.detect$object), 
                 "data"         = data,
                 "region.table" = region.table, 
                 "sample.table" = sample.table)
  
  ## survey result plots
  if (plot == TRUE) {
    plot(survey.results)
  }
  
  return(tables)
  
}