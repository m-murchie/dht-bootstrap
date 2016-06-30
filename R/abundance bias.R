source("R/double.obs.sim.mod.R")
source("R/boot.dht.R")

detect.1 <- make.detectability(key.function = "hn", scale.param = 20,
                               truncation = 20)

detect.2 <- make.detectability(key.function = "hn", scale.param = 13.5,
                               truncation = 20)

ddf.analyses <- make.ddf.analysis.list(dsmodel = list(~cds(key = "hn",formula = 
                                                             ~distance*observer)),
                                       method = "ds", truncation = 20)

pop.description <- make.population.description(N = 500, density.obj = pop.density,
                                               region = region, fixed.N = TRUE)

tables <- double.obs.sim.mod(region, survey.design, pop.description, detect.1,
                             detect.2, ddf.analyses)


ddf.model <- ddf(method = 'io.fi', mrmodel=~glm(link='logit', formula=~distance*
                                                                       observer), 
                 data = tables$data, meta.data=list(width=20))
dht.results <- dht(ddf.model, tables$region.table, tables$sample.table, 
                   tables$obs.table)


boot.dht(tables, B=49, trunc=20, hist = TRUE)

