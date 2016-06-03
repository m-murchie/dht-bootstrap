#           --- CRABEATER SEAL DATA EXAMPLE: BOOTSTRAP VARIANCE ESTIMATION ---

## load mrds library and crabeaster seal data
rm(list = ls())
library(mrds)
library(Distance)
crabseal <- read.csv("crabbieMRDS.csv")
crabseal$Sample.Label <- as.numeric(crabseal$Sample.Label)   # number transects

crab.ddf.io <- ddf(method="trial.fi",
                   mrmodel=~glm(link="logit", formula=~distance),
                   data=crabseal, meta.data=list(width=700))
summary(crab.ddf.io)["Nhat"]$Nhat


## bootstrap
B <- 999
strip.count <- max(crabseal$Sample.Label)   # number of transects 
N.hats <- NULL   

for (i in 1:B) {
  strip.sample <- sample(1:strip.count, strip.count, replace=TRUE)
  boot.sample <- crabseal[crabseal$Sample.Label==strip.sample[1],]
  for (j in 2:strip.count) {
    boot.sample <- rbind(boot.sample, crabseal[crabseal$Sample.Label==strip.sample[j],])
  }
  crab.ddf <- ddf(method="trial.fi",
                  mrmodel=~glm(link="logit", formula=~distance),
                  data=boot.sample, meta.data=list(width=700))
  N.hats[i] <- summary(crab.ddf)[9]$Nhat   # store abundance estimate
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

