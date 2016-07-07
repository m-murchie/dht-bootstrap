n <- 499
N.c <- vector(length=n)
N.hats <- vector(length=n)
diff <- vector(length=n)


for (i in 1:n) {
  lines <- create.transects(n=9, trunc = 20)
  
  unlink("shapes/*")
  write.transects(lines, "shapes")
  survey.design <- make.design(transect.type = "Line",
                               design.details = c("Parallel","Systematic"), region = region,
                               plus.sampling =FALSE, path = shapefile.pathway)
  
  tables <- double.obs.sim.mod.2(region, survey.design, pop.description, detect.1,
                                 detect.2, ddf.analyses, seed=runif(1,max=9999999))
  
  
  ddf.model <- ddf(method = 'io.fi', mrmodel=~glm(link='logit', formula=~distance*
                                                    observer), 
                   data = tables$data, meta.data=list(width=20))
  dht.results <- dht(ddf.model, tables$region.table, tables$sample.table, 
                     subset=1==1)
  
  N.c[i] <- ddf.model$Nhat
  N.hats[i] <- dht.results$individuals$N$Estimate
  diff[i] <- tables$true.Nc - ddf.model$Nhat
}N
