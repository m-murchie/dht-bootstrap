source("R/double.obs.sim.R")
source("R/bootstrap function.R")


detect.1 <- make.detectability(key.function = "hn", scale.param = 50,
                               truncation = 30)

detect.2 <- make.detectability(key.function = "hn", scale.param = 40,
                               truncation = 30)

tables <- double.obs.sim(region, survey.design, pop.description, detect.1, detect.2, 
                         ddf.analyses, plot=TRUE)

tables$data["Region.Label"] <- rep("study.area", length(tables$data$object))


ddf.model <- ddf(method = 'io.fi', mrmodel=~glm(link='logit', formula=~distance), 
                 data = tables$data, meta.data=list(width=30))
dht.results <- dht(ddf.model, tables$region.table, tables$sample.table, 
                   tables$obs.table)


boot.dht(tables, trunc=30, hist=TRUE)
