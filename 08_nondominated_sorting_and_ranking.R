require("dplyr")
require("emoa")
require("mco")
require("ggplot2")
require("rPref")


file_name <- paste('./data/predictions_smith_waterman_PAM30_2021-06-08_14_22_03.csv', sep="")
#file_name <- paste('./data/predictions_smith_waterman_BLOSUM45_2021-06-08_15_57_07.csv', sep="")
df <- read.csv(file = file_name)

hist(df$pair_align_score)

df_objectives <- df %>%
  mutate(probability_inv = 1 - Probability,
         distance_neg = -pair_align_score) %>%
  select(probability_inv, distance_neg) #%>% head() 

View(df)


# First pareto front
####################

m_objectives <- t(as.matrix(df_objectives))

nondominated <- nondominated_points(m_objectives)

rank <- nds_rank(nondominated) #, partial)

dim(m_objectives)
dim(nondominated)
dim(rank)

nondominated
rank

# plots Pareto fronts 
m_objectives_t <- as.data.frame(t(nondominated))
plot(m_objectives_t$probability_inv, m_objectives_t$distance_neg, 
     xlim = c(min(df_objectives$probability_inv), max(df_objectives$probability_inv)),
     ylim = c(min(df_objectives$distance_neg), max(df_objectives$distance_neg)))
points(m_objectives_t$probability_inv, m_objectives_t$distance_neg, col="blue", pch="o")
lines(m_objectives_t$probability_inv, m_objectives_t$distance_neg, col="blue", lty=1)


# Additional pareto fronts
##########################

for (i in 1:20) {
  # https://www.datasciencemadesimple.com/join-in-r-merge-in-r/
  df_objectives <- df_objectives %>% 
    anti_join(as.data.frame(t(nondominated)), by="probability_inv")#, distance_neg")

  m_objectives <- t(as.matrix(df_objectives))
  
  nondominated <- nondominated_points(m_objectives)
  
  rank <- nds_rank(nondominated) #, partial)
  
  dim(m_objectives)
  dim(nondominated)
  dim(rank)
  
  nondominated
  rank
  
  color <- "green"
  if (!i %% 2) {color <- "blue"}
  
  # plots Pareto fronts 
  m_objectives_t <- as.data.frame(t(nondominated))
  points(m_objectives_t$probability_inv, m_objectives_t$distance_neg, col=color, pch="o")
  lines(m_objectives_t$probability_inv, m_objectives_t$distance_neg, col=color, lty=1)
}  

