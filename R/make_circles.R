# functions for creating circular polygons of specified size
make_circle <- function(
	n = 100, # number of points, equally spaced around the edge of a circle
	k,  # scaling value (radius) - larger for bigger circles. defaults to 1
  offset = c(0, 0)
){
	if(missing(k)){k <- 1}
	# get coordinates
	theta <- (2 * (pi/(n-1)) * seq_len((n-1)) ) # -alpha
	values <- data.frame(
		x = k * cos(theta),
		y = k * sin(theta))
	rownames(values) <- NULL
  values$x <- values$x + offset[[1]]
  values$y <- values$y + offset[[2]]
  # return(values)
	return(values[c(seq_len(n-1), 1), ])
}
# note the result has nrow == n, but n-1 unique rows

# run make circle in a loop
make_circles <- function(
  df, # should be returned by packcircles::circleProgressiveLayout
  gap = 0
){
  df_list <- split(df, seq_len(nrow(df)))
  result <- as.data.frame(do.call(rbind,
    lapply(df_list, function(a){
      make_circle(k = a$radius - gap, offset = c(a$x, a$y))
    })
  ))
  result$group <- rep(seq_along(df_list), each = 100)
  return(result)
}