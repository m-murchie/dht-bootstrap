create.transects <- function(n, trunc) {
  
  n_segs <- 2
  space <- 1000/n
  start <- runif(1, min = trunc, max = (space - trunc))
  
  y.coords <- seq(start, 1000, space)
  legs <- as.character(1:n)
  
  lines <- data.frame(x   = c(rep(seq(0, 1000, len = n_segs), n)),
                      y   = c(rep(y.coords, each = n_segs)),
                      leg = c(rep(legs, each = n_segs)))
  
  return(lines)

}
