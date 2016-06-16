source("R/double.obs.sim.R")

detect.1 <- make.detectability(key.function = "hn", scale.param = 80,
                               truncation = 30)

detect.2 <- make.detectability(key.function = "hn", scale.param = 60,
                               truncation = 30)

tables <- double.obs.sim(region, survey.design, pop.description, detect.1, detect.2, 
                         ddf.analyses, plot=TRUE, seed=999)


ddf.model <- ddf(method = 'io.fi', mrmodel=~glm(link='logit', formula=~distance), 
                 data = tables$data, meta.data=list(width=30))
dht.results <- dht(ddf.model, tables$region.table, tables$sample.table, 
                   tables$obs.table)

B <- 999
lines <- unique(tables$obs.table$Sample.Label)
lines.count <- length(lines)   # number of transects
width = 30
e.rates <- NULL   # encounter rates  
f.zero <- NULL   # f(0), pdf of observed distance evaluated at 0 distance
N.hats <- NULL   # abundance estimates

for (i in 1:B) {
  line.sample <- sample(lines, lines.count, replace=TRUE)
  boot.data <- tables$data[tables$data$Sample.Label==line.sample[1],]
  boot.sample <- tables$sample.table[tables$sample.table$Sample.Label==
                                       line.sample[1],]
  for (j in 2:lines.count) {
    boot.data <- rbind(boot.data, tables$data[tables$data$Sample.Label
                                                  ==line.sample[j],])
    boot.sample <- rbind(boot.sample, tables$sample.table
                               [tables$sample.table$Sample.Label==line.sample[j],])
  }
  boot.ddf <- ddf(method = 'io.fi', mrmodel=~glm(link='logit', formula=~distance), 
                   data = boot.data, meta.data=list(width=30))
  boot.effort <- sum(boot.sample$Effort)
  covered.area <- 2*boot.effort*width
  region.area <- tables$region.table$Area
  N.hats[i] <- boot.ddf$Nhat*region.area/covered.area
  e.rates[i] <- summary(boot.ddf)$n/boot.effort
}


boot.mean.er <- mean(e.rates)   # mean of bootstrap encounter rate estimates
boot.se.er <- sd(e.rates)   # estimate standard error from bootstrap estimates
boot.CV.er <- boot.se.er/boot.mean.er
c(boot.mean.er, boot.se.er, boot.CV.er)  # encounter rate, se and CV estimates

hist(e.rates)  # histogram of bootstrap encounter rate estimates
hist(N.hats)   # histogram of bootstrap abundance estimates
