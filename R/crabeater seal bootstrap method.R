#           --- CRABEATER SEAL DATA EXAMPLE: BOOTSTRAP VARIANCE ESTIMATION ---

## load mrds library and crabeaster seal data
rm(list = ls())
library(mrds)
crabseal <- read.csv("crabbieMRDS.csv")
crabseal$Sample.Label <- as.numeric(crabseal$Sample.Label)   # number transects

## detection function
crab.ddf <- ddf(method="trial.fi",
                   mrmodel=~glm(link="logit", formula=~distance),
                   data=crabseal, meta.data=list(width=700))

## estimate abundance and variance
dht(region=tables$region.table, sample=tables$sample.table, obs=tables$obs.table,
    model=crab.ddf, options=list(convert.units=0.001))

## bootstrap
B <- 100
strip.count <- max(crabseal$Sample.Label)   # number of transects 
N.hats <- NULL   

for (i in 1:B) {
  strip.sample <- sample(1:strip.count, strip.count, replace=TRUE)
  boot.sample <- crabseal[crabseal$Sample.Label==strip.sample[1],]
  dht.sample <- tables$sample.table[tables$sample.table$Sample.Label %in% strip.sample,]
  for (j in 2:strip.count) {
    boot.sample <- rbind(boot.sample, crabseal[crabseal$Sample.Label==strip.sample[j],])
  }
  boot.ddf <- ddf(method="trial.fi",
                  mrmodel=~glm(link="logit", formula=~distance),
                  data=boot.sample, meta.data=list(width=700))
  boot.dht <- dht(region=tables$region.table, sample=dht.sample, obs=tables$obs.table,
                  model=boot.ddf, options=list(convert.units=0.001))
  N.hats[i] <- boot.dht$cluster$N[1,2]   # store abundance estimate
}

# abundance of original data
crab.ddf.io <- ddf(method="trial.fi",
                   mrmodel=~glm(link="logit", formula=~distance),
                   data=crabseal, meta.data=list(width=700))
N.hats[1000] <- summary(crab.ddf.io)[9]$Nhat

boot.mean.N <- mean(N.hats)   # mean of bootstrap abundance estimates
boot.se.N <- sd(N.hats)   # estimate standard error from bootstrap estimates
boot.CV.N <- sd(N.hats)/mean(N.hats)
c(boot.mean.N, boot.se.N, boot.CV.N)  # adundance, se and CV estimates

hist(N.hats)  # histogram of bootstrap abundance estimates
#---

