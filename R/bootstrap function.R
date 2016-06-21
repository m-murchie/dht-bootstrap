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
  
  
  if (!("Sample.Label" %in% names(tables$data)) | 
      !("Region.Label" %in% names(tables$data))) {
        stop("Data must contain Sample.Label and Region.Label column 
             to execute bootstrap")
  }
  
  
  # bootstrap estimates of abundance
  lines <- unique(tables$obs.table$Sample.Label)
  lines.count <- length(lines)   # number of transects
  e.rates <- NULL   # encounter rates  
  N.hats <- NULL   # abundance estimates
  
  
  for (i in 1:B) {
    line.sample <- sample(lines, lines.count, replace=TRUE)
    boot.data <- tables$data[tables$data$Sample.Label==line.sample[1],]
    boot.data$Sample.Label <- rep(1, length(boot.data$Sample.Label))
    boot.sample.table <- tables$sample.table[tables$sample.table$Sample.Label
                                             ==line.sample[1],]
  
    for (j in 2:lines.count) {
      line.data <- tables$data[tables$data$Sample.Label==line.sample[j],]
      line.data$Sample.Label <- rep(j, length(line.data$Sample.Label))
      boot.data <- rbind(boot.data, line.data)
      boot.sample.table <- rbind(boot.sample.table, 
                                 tables$sample.table[tables$sample.table$
                                                     Sample.Label==line.sample[j],])
    }
  
    boot.objects <- c(1:(length(boot.data$object)*0.5))
    boot.objects <- rep(boot.objects, each=2)
    boot.data$object <- boot.objects
    boot.sample.table$Sample.Label <- c(1:lines.count)
    
    boot.ddf <- ddf(method = 'io.fi', mrmodel=~glm(link='logit', formula=~distance), 
                    data = boot.data, meta.data=list(width=trunc))
    boot.dht <- dht(boot.ddf, tables$region.table, boot.sample.table, 
                    subset=1==1)
    
    boot.effort <- boot.dht$individuals$summary$Effort
    N.hats[i] <- boot.dht$individuals$N$Estimate
    e.rates[i] <- summary(boot.ddf)$n/boot.effort
  }
  
  
  # original data
  original.ddf <- ddf(method = 'io.fi', mrmodel=~glm(link='logit', 
                                                     formula=~distance),
                      data = tables$data, meta.data=list(width=trunc))
  original.dht <- dht(original.ddf, tables$region.table, tables$sample.table, 
                     tables$obs.table)
  N.hats[B+1] <- original.dht$individuals$N$Estimate
  
  
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
  boot.result <- data.frame(Method   = "Bootstrap",
                            Label    = "Total",
                            Estimate = boot.mean,
                            se       = boot.se,
                            cv       = boot.cv,
                            lcl      = N.hats[(B+1)*0.025],
                            ucl      = N.hats[(B+1)*0.975])
  
  original.result <- original.dht$individuals$N
  original.result <- cbind(Method = "Original", original.result[1:6])
  
  result <- list(original  = original.result,
                 bootstrap = boot.result)

  return(result)

}


#---
