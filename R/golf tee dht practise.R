#           --- GOLF DATA EXAMPLE: BOOTSTRAP VARIANCE ESTIMATION ---

## load mrds library and golf tee data
rm(list = ls())
library(mrds)
data("book.tee.data")

## tidy up tee data
tee.data <- book.tee.data$book.tee.dataframe
tee.data$sex <- as.factor(tee.data$sex)             ## sex and exposure are factor
tee.data$exposure <- as.factor(tee.data$exposure)   ## variables
tee.region <- book.tee.data$book.tee.region
tee.samples <- book.tee.data$book.tee.samples
tee.obs <- book.tee.data$book.tee.obs

## add sample label column to tee data
tee.obs <- tee.obs[order(tee.obs$object),]
Sample.Label <- tee.obs$Sample.Label
Sample.Label <- rep(Sample.Label, each=2)
tee.data["Sample.Label"] <- Sample.Label

## add region label column to tee data
Region.Label <- tee.obs$Region.Label
Region.Label <- rep(Region.Label, each=2)
tee.data["Region.Label"] <- Region.Label

## bootstrap sample for tee data
B <- 999
region.strips <- unique(tee.data[tee.data$Region.Label==2,]$Sample.Label)
strip.count <- length(region.strips)   # number of transects 
N.hats <- NULL                                                    

for (i in 1:B) {
  strip.sample <- sample(region.strips, strip.count, replace=TRUE)
  boot.sample <- tee.data[tee.data$Sample.Label==strip.sample[1],]
  for (j in 2:strip.count) {
    boot.sample <- rbind(boot.sample, tee.data[tee.data$Sample.Label==strip.sample[j],])
  }
  tee.model <- ddf(method = 'trial.fi', mrmodel=~glm(link='logit', formula=~distance+exposure), 
                   data = boot.sample, meta.data=list(width=4))
  N.hats[i] <- summary(tee.model)["Nhat"]$Nhat
}

N## abundance estimate from original sample
tee.model.orig <- ddf(method = 'trial.fi', mrmodel=~glm(link='logit', formula=~distance+exposure), 
                      data = tee.data, meta.data=list(width=4))
N.hats[1000] <- summary(tee.model.orig)["Nhat"]$Nhat

boot.mean.N <- mean(N.hats)   # mean of bootstrap abundance estimates
boot.se.N <- sd(N.hats)   # estimate standard error from bootstrap estimates
boot.CV.N <- sd(N.hats)/mean(N.hats)
c(boot.mean.N, boot.se.N, boot.CV.N)  # adundance, se and CV estimates

hist(N.hats)  # histogram of bootstrap abundance estimates
#---

## NOTE: se and CV estimates obtained from the bootstrap are larger than
##       those calculated using the delta method - tends to underestimate se.