create.transects <- function(n) {
  
  n_segs <- 2
  
  y.coords <- seq(1000/(n+1), 1000*n/(n+1), 1000/(n+1))
  legs <- as.character(1:n)
  
  lines <- data.frame(x   = c(rep(seq(0, 1000, len = n_segs), n)),
                      y   = c(rep(y.coords, each = n_segs)),
                      leg = c(rep(legs, each = n_segs)))
  
  return(lines)

}
