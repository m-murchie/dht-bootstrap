#           --- GOLF DATA EXAMPLE: DHT VARIANCE ESTIMATION ---

## load mrds library and golf tee data
rm(list = ls())
library(mrds)
data("book.tee.data")

## tidy up tee data
tee.data <- book.tee.data$book.tee.dataframe
tee.data$sex <- as.factor(tee.data$sex)        ## sex and exposure are factor
tee.data$exposure <- as.factor(tee.data$exposure)   ## variables
tee.region <- book.tee.data$book.tee.region
tee.samples <- book.tee.data$book.tee.samples
tee.obs <- book.tee.data$book.tee.obs

## fit dection function
ddf.model <- ddf(method = 'trial.fi', mrmodel=~glm(link='logit', formula=~distance+size+exposure), 
                data = tee.data, meta.data=list(width=4))

summary(ddf.model)["Nhat"]$Nhat

## GOF and qq-plot
ddf.gof(ddf.model, main="Detection Function QQ-plot")

## create table for data collected by two observer approach
tee.det.table <- det.tables(ddf.model)
plot(tee.det.table) 


## density and variance estimation 
dht(ddf.model, tee.region, tee.samples, tee.obs)

#---