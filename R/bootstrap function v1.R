#                         --- BOOTSTRAP FUNCTION ---


boot.dht <- function(tables, B = 999, trunc, hist = FALSE) {
# INPUTS:    data    dataframe containing data to be analysed, WITH Sample.Label
#                    column of numeric values 
#            B       number of bootstrap replicates
#            hist    produce histogram of abundance estimates and encounter rates
#            trunc   truncation distance, necessary for meta.data 
# OUTPUTS:   summary table for bootstrap abundance estimates
#            histogram of abundance estimates and encounter rates, if hist = TRUE
# FUNCTIONS: ddf
  if (!("Sample.Label" %in% names(tables$data))) {
    stop("Data must contain Sample.Label column to execute bootstrap")
  }
  
  # bootstrap estimates of abundance
  lines <- unique(tables$obs.table$Sample.Label)
  lines.count <- length(lines)   # number of transects
  e.rates <- NULL   # encounter rates  
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
                    data = boot.data, meta.data=list(width=trunc))
    boot.effort <- sum(boot.sample$Effort)
    covered.area <- 2*boot.effort*trunc
    region.area <- tables$region.table$Area
    N.hats[i] <- boot.ddf$Nhat*region.area/covered.area
    e.rates[i] <- summary(boot.ddf)$n/boot.effort
  }
  
  
  # original data
  original.effort <- sum(tables$sample.table$Effort)
  original.ddf <- ddf(method = 'io.fi', mrmodel=~glm(link='logit', 
                                                     formula=~distance), 
                      data = boot.data, meta.data=list(width=trunc))
  dht.results <- dht(original.ddf, tables$region.table, tables$sample.table, 
                     tables$obs.table)
  N.hats[B+1] <- dht.results$individuals$N$Estimate[[1]]
  e.rates[B+1] <- summary(original.ddf)$n/original.effort
  
  
  # summary characteristics
  N.hats <- sort(N.hats)
  boot.mean <- mean(N.hats)
  boot.se <- sd(N.hats)   
  boot.cv <- boot.se/boot.mean  
  
  
  # histogram of bootstap abundance estimates and encounter rates
  if (hist == TRUE) {
    par(mfrow=c(2,2))
    hist(N.hats)
    hist(e.rates)
    par(mfrow=c(1,1))
  }
  
  
  # results
  result <- data.frame(Estimate = boot.mean,
                       SE       = boot.se,
                       CV       = boot.cv,
                       lcl      = N.hats[(B+1)*0.025],
                       ucl      = N.hats[(B+1)*0.975])
  rownames(result) <- c("Bootstrap")
  
  return(result)
}


#---
