## taken from DL Miller, spatlaugh repo

# write out a transect as a shapefile

# this is a completely Byzantine process

library(sp)
library(rgdal)

# transect should have the following columns:
#  x,y   x,y coordinates of **endpoints**
#  leg   survey leg associated
# file   filename (well actually foldername) to write to
write.transects <- function(transect, file){
  
  llist <- list()
  
  Sample.Label <- c()
  
  # make each leg a different Line object
  for(this_leg in unique(transect$leg)){
    
    # get the transect bits for this leg
    tr <- transect[transect$leg==this_leg,]
    
    
    for(i in 1:(nrow(tr)-1)){
      this_label <- paste0(this_leg, "-", i)
      ll <- Line(tr[,c("x","y")][i:(i+1),])
      llist <- c(llist, Lines(ll, ID=this_label))
      Sample.Label <- c(Sample.Label, this_label)
    }
  }
  
  
  ll <- SpatialLines(llist)
  
  dat <- data.frame(Sample.Label=Sample.Label)
  rownames(dat) <- Sample.Label
  
  ll <- SpatialLinesDataFrame(ll, data=dat)
  
  writeOGR(ll, file, "data", "ESRI Shapefile" )
}