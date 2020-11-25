# build data for plotting taxonomic information
# solution here is based on circle packing

taxon_plot_data <- function(
  labels, # vector of taxon labels
  counts, # vector of counts
  color_scale = "viridis",
  log_scale = TRUE, # logical - should the result be log + 1 -scaled?
  gap = 0 # how much should the radius be reduced to decrease circle overlap?
){
  if(log_scale){
    taxon_counts <- log(counts + 1)
  }else{
    taxon_counts <- counts
  }
  circle_df <- packcircles::circleProgressiveLayout(taxon_counts)
  circle_plot_df <- make_circles(circle_df, gap = gap)
  circle_plot_df$label <- factor(circle_plot_df$group,
    levels = unique(circle_plot_df$group),
    labels = labels)
  circle_plot_df$count <- rep(counts, each = 100)
  circle_plot_df$plotly_text <- paste0(
    "<b>",
    circle_plot_df$label,
    "</b><br>n = ",
    formatC(circle_plot_df$count, big.mark = ","))
  circle_plot_df$color <- rep(
    do.call(color_scale, args = list(
      n = length(labels),
      direction = -1)),
    each = 100)
  return(circle_plot_df)
}