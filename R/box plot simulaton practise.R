#                         --- BOX PLOT SIMULATION ---

# load parallel library
library(parallel)

n <- 10

seeds <- runif(n)*1e7

# Calculate the number of cores
no.cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no.cores)
clusterEvalQ(cl, library(DSsim))
clusterExport(cl, c('region', 'survey.design', 'pop.description', 
                    'detect.1', 'detect.2', 'ddf.analyses', 'double.obs.sim.mod'))


# create tables
store.tables <- list(length = n)
store.tables <- parLapply(cl, seeds, function(x) 
                          double.obs.sim.mod(region, survey.design, pop.description, 
                                             detect.1, detect.2, ddf.analyses, 
                                             seed = x))

stopCluster(cl)

#---

# Initiate cluster
cl <- makeCluster(no.cores)
clusterEvalQ(cl, library(mrds))
clusterExport(cl, 'boot.dht')

results <- list(length = n)
results <- parLapply(cl, store.tables, function(x)
                                       boot.dht(x, B=49, trunc=20))


stopCluster(cl)

#---

library(plotly)
orig.se <- unlist(lapply(results, function(x) x$original$se))
boot.se <- unlist(lapply(results, function(x) x$bootstrap$se))

plot_ly(y = orig.se, type = "box", name="original") %>%
        add_trace(y = boot.se, type = "box", name="bootstrap") %>% 
        layout(xaxis = list(title = "method"), yaxis = list(title = "se"))
