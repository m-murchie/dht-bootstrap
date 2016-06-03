#                    --- BOOTSTRAP FUNCTION: VERSION 1 ---

rm(list = ls())

boot <- function(data, B = 999, hist = FALSE, trunc) {
# INPUTS:    data    dataframe containing data to be analysed, WITH Sample.Label
#                    column of numeric values 
#            B       number of bootstrap replicates
#            hist    produce histogram of abundance estimates
#            trunc   truncation distance, necessary for meta.data 
# OUTPUTS:   tables comparing bootstrap and original sample estimates
#            histogram of abundance estimates, if hist = TRUE
# FUNCTIONS: ddf
  if (!("Sample.Label" %in% names(data))) {
    stop("Data must contain Sample.Label column to execute bootstrap")
  }
  
  # bootstrap estimates of abundance
  strip.count <- max(data$Sample.Label)   # number of transects 
  N.hats <- NULL 
  
  for (i in 1:B) {
    strip.sample <- sample(1:strip.count, strip.count, replace=TRUE)
    boot.sample <- data[data$Sample.Label==strip.sample[1],]
    for (j in 2:strip.count) {
      boot.sample <- rbind(boot.sample, data[data$Sample.Label==strip.sample[j],])
    }
    ddf.model <- ddf(method="trial.fi",
                     mrmodel=~glm(link="logit", formula=~distance),
                     data=boot.sample, meta.data = list(width=trunc))
    N.hats[i] <- summary(ddf.model)["Nhat"]$Nhat   # store abundance estimate
  }
  
  # abundance estimate of original sample 
  ddf.model.orig <- ddf(method="trial.fi",
                        mrmodel=~glm(link="logit", formula=~distance),
                        data=data, meta.data = list(width=trunc))
  N.hats[B+1] <- summary(ddf.model.orig)["Nhat"]$Nhat
  
  # histogram of abundance estimates
  if (hist == TRUE) {
    hist(N.hats)
  }
  
  # summary characteristics
  boot.mean <- mean(N.hats)
  boot.se <- sd(N.hats)   
  boot.CV <- boot.se/boot.mean  
  
  sample.mean <- summary(ddf.model.orig)["Nhat"]$Nhat
  sample.se <- summary(ddf.model.orig)["Nhat.se"]$Nhat.se
  sample.CV <- sample.se/sample.mean

  result.mean <- c(boot.mean, sample.mean)
  result.se <- c(boot.se, sample.se)
  result.CV <- c(boot.CV, sample.CV)
  
  result <- data.frame(Estimate = result.mean,
                       SE       = result.se,
                       CV       = result.CV)
  rownames(result) <- c("Bootstrap", "Original")
  
  return(result)
}

#---

## GOLF TEE EXAMPLE
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

boot(tee.data, hist=TRUE, trunc = 4)

#---

## CRABEATER SEAL EXAMPLE

library(Distance)
crabseal <- read.csv("crabbieMRDS.csv")
crabseal$Sample.Label <- as.numeric(crabseal$Sample.Label)   # number transects

boot(crabseal, hist=TRUE, trunc = 700)

