require("dplyr")
require("emoa")
require("mco")
require("ggplot2")
library("plotly")

png('./output/pareto_fronts_PAM30.png', width=2000, height=1600, res=300)

file_name <- './data/predictions_smith_waterman_PAM30_2021-06-16_23_18_07.csv'
                
df <- read.csv(file=file_name)

df_complemented_negated <- df %>%
  mutate(probability_complement = 1 - Probability,
         distance_neg = -pair_align_score) %>%
  select(Protein, Sequence, probability_complement, distance_neg, Info_peptide, Info_PepID) #%>% head() 

df_objectives <- df_complemented_negated %>%
  select(probability_complement, distance_neg) #%>% head() 

m_objectives <- t(as.matrix(df_objectives))

nondominated <- nondominated_points(m_objectives)

rank <- nds_rank(m_objectives) #, partial)

file_name_split <- unlist(strsplit(file_name, '_'))

# Prepare plot for Pareto fronts 
plot(c(), c(), 
     xlim=c(min(df_objectives$probability_complement), max(df_objectives$probability_complement)),
     ylim=c(min(df_objectives$distance_neg), max(df_objectives$distance_neg)),
     main=file_name_split[4],
     xlab='Probability Complement',
     ylab='Negated Pairwise Alignment Score')

# Plot Pareto fronts
for (idx in 1:max(rank)) {
  color <- 'blue'
  shape <- 1
  if (!idx %% 2) {
    color <- 'green'
    shape <- 4
  }
  
  # Prepare Pareto points
  pareto_front_points_indicies <- which(rank %in% c(idx))
  df_pareto_front_points <- df_objectives[pareto_front_points_indicies, 
                                          c('probability_complement', 'distance_neg')]
    
  # plots Pareto fronts 
  points(df_pareto_front_points$probability_complement, 
         df_pareto_front_points$distance_neg, 
         pch=shape)
  lines(df_pareto_front_points$probability_complement, 
        df_pareto_front_points$distance_neg, 
        col=color, lty=1)
}

dev.off()



# ######################################################################################
# 
# # https://stat.ethz.ch/R-manual/R-devel/library/graphics/html/identify.html
# 
# ## A function to use identify to select points, and overplot the
# ## points with another symbol as they are selected
# identifyPch <- function(x, y = NULL, n = length(x), plot = FALSE, pch = 19, ...)
# {
#   xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
#   sel <- rep(FALSE, length(x))
#   while(sum(sel) < n) {
#     ans <- identify(x[!sel], y[!sel], labels = which(!sel), n = 1, plot = plot, ...)
#     if(!length(ans)) break
#     ans <- which(!sel)[ans]
#     points(x[ans], y[ans], pch = pch)
#     sel[ans] <- TRUE
#   }
#   ## return indices of selected points
#   which(sel)
# }
# 
# if(dev.interactive()) { ## use it
#   x <- rnorm(50); y <- rnorm(50)
#   plot(x,y); identifyPch(x,y) # how fast to get all?
# }


####################################################################################

# https://cengel.github.io/R-data-viz/interactive-graphs.html
# https://stackoverflow.com/a/34617687

#plot_ly(df_pareto_front_points, x = ~probability_complement, y = ~distance_neg)

p <- ggplot(data=df_objectives, aes(x=probability_complement, y=distance_neg)) +
  geom_point(aes(text=sprintf("Protein ID: %s Sequence: %s<br>Info_PepID: %s Info_peptide: %s", 
                              Protein, Sequence, 
                              Info_PepID, Info_peptide)))

ggplotly(p) 

####################################################################################



