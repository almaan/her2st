require(ggplot2)
require(RColorBrewer)

# FUNCTIONS -----------

make_empty_theme <- function(){
  theme_empty <- theme_bw()
  theme_empty$line <- element_blank()
  theme_empty$rect <- element_blank()
  theme_empty$strip.text <- element_blank()
  theme_empty$axis.text <- element_blank()
  theme_empty$plot.title <- element_text(colour = "white",
                                         family = "calibri",
                                         face = "bold",
                                         size = 20,
                                         hjust = 0,
                                         margin = margin(t = 10, b = -20)) #element_blank()
  theme_empty$axis.title <- element_blank()
  theme_empty$legend.position <- c(0.95, 0.5)
  theme_empty$legend.background <- element_rect(fill = "transparent",
                                                colour = "transparent")
  theme_empty$legend.text <- element_text(colour = "white")
  theme_empty$legend.title <- element_text(colour = "white")
  theme_empty$strip.background <- element_rect(fill="black")
  theme_empty$legend.key.size <- unit(10,"mm")
  theme_empty$legend.margin <- margin(10,50,10,50)

  return(theme_empty)

}


## VARIABLES ------------


COLORS <- list("Green" = RColorBrewer::brewer.pal(n = 9, name = "Greens"),
               "Blue" = RColorBrewer::brewer.pal(n = 9, name = "Blues"),
               "Red" = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
               "Blue-Yellow-Red" = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")),
               "Spectral" = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")),
               "viridis" = viridis::viridis(n = 9),
               "magma" = viridis::magma(n = 9)
               )

EDGESTROKES <- list("On" = 0.5,
                   "Off" = 0.0)


