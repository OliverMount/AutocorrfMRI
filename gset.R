gset <- function(g)
{
  # Edge and vertex Setting for the graph
  E(g)$color = "black"
  E(g)$width = 2
  E(g)$arrow.width = .25
  E(g)$label.cex = 1.4
  E(g)$label.color = "black"
  V(g)$label.color = "black"
  V(g)$size =30
  return(g)
}