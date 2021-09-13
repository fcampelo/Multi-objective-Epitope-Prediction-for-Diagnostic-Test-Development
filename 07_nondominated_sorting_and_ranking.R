require("dplyr")
require("emoa")
require("mco")
require("ggplot2")
require("rPref")


#file_name <- paste('./data/predictions_smith_waterman_PAM30_2021-06-08_14_22_03.csv', sep="")
file_name <- paste('./data/predictions_smith_waterman_BLOSUM45_2021-06-08_15_57_07.csv', sep="")
df <- read.csv(file = file_name)

hist(df$pair_align_score)

df_objectives <- df %>%
  mutate(probability_inv = 1 - Probability,
         distance_neg = -pair_align_score) %>%
  select(probability_inv, distance_neg) #%>% head() 

#print(df_objectives)  
#print(t(as.matrix(df_objectives)))

m_objectives <- t(as.matrix(df_objectives))

nondominated <- nondominated_points(m_objectives)

rank <- nds_rank(nondominated) #, partial)

dim(m_objectives)
dim(nondominated)
dim(rank)

nondominated
rank

# #plt <- qplot(hypervolume, data=subset(df, hypervolume > 127), fill=method, geom="density")
# plt <- qplot(nondominated, fill=method, geom="density")
# plt



# plots Pareto fronts 
show_front <- function(pref) {
  #m_objectives_t <- as.data.frame(t(m_objectives))
  m_objectives_t <- as.data.frame(t(nondominated))
  #plot(m_objectives_t$probability_inv, m_objectives_t$distance_neg)
  #sky <- psel(m_objectives_t, pref)
  #plot_front(m_objectives_t, pref, col = rgb(0, 0, 1))
  #points(sky$probability_inv, sky$distance_neg, lwd = 1)
  points(m_objectives_t$probability_inv, m_objectives_t$distance_neg, col="blue", pch="o")
  lines(m_objectives_t$probability_inv, m_objectives_t$distance_neg, col="blue",lty=1)
  
}

# do this for all four combinations of Pareto compositions
show_front(low(probability_inv)  * low(distance_neg))


################################################################################

# df_objectives <- df_objectives %>%
#   filter(probability_inv!=0.135166 & distance_neg!=-62.0)  %>%
#   filter(probability_inv!=0.135166 & distance_neg!=-62.0)  %>%
#   filter(probability_inv!=0.135166 & distance_neg!=-62.0)  %>%
#   filter(probability_inv!=0.1728069 & distance_neg!=-91.0)  %>%
#   filter(probability_inv!=0.1986627 & distance_neg!=-170.0)  %>%
#   filter(probability_inv!=0.2766 & distance_neg!=-187.0)
# 
df_objectives <- df_objectives %>% 
  anti_join(nondominated, by="probability_inv, distance_neg")

  
m_objectives <- t(as.matrix(df_objectives))

nondominated <- nondominated_points(m_objectives)

rank <- nds_rank(nondominated) #, partial)

dim(m_objectives)
dim(nondominated)
dim(rank)

nondominated
rank


# plots Pareto fronts 
show_front_2 <- function(pref) {
  #m_objectives_t <- as.data.frame(t(m_objectives))
  m_objectives_t <- as.data.frame(t(nondominated))
  points(m_objectives_t$probability_inv, m_objectives_t$distance_neg, col="green", pch="o")
  lines(m_objectives_t$probability_inv, m_objectives_t$distance_neg, col="green",lty=1)
}

# do this for all four combinations of Pareto compositions
show_front_2(low(probability_inv)  * low(distance_neg))




#  |

# https://cran.r-project.org/web/packages/emoa/emoa.pdf
# http://www.statistik.tu-dortmund.de/~olafm/software/emoa/
# https://p-value.net/
# https://github.com/olafmersmann/emoa

# nondominated_points Nondominated points
#   Description
#     Return those points which are not dominated by another point in points. This is the Pareto front
#     approximation of the point set.
#   Usage
#     nondominated_points(points)
#   Arguments
#     points - Matrix of points, one per column.
#   Value
#     Those points in points which are not dominated by another point.

# nds_rank Nondominated sorting ranks
#   Description
#     Perform (partial) nondominated sort of the points in points and return the rank of each point.
#   Usage
#     nds_rank(points, partial)
#     nondominated_ordering(points, partial)
#   Arguments
#     points - Matrix containing points one per column.
#     partial - Optional integer specifying the number of points for which the rank should be
#               calculated. Defaults to all points.
#   Value
#     Vector containing the ranks of the first partial individuals or all individuals.

# e.g. nondominated_points
# https://github.com/olafmersmann/emoa/blob/master/examples/ex-benchmark.R
#
# e.g. nds_rank
# https://github.com/olafmersmann/emoa/search?q=nds_rank




#source("sms_emoa.r")
#
# sms_emoa <- function(f, lower, upper, ...,
#                      control=list(mu=100L,
#                                   sbx.n=15, sbx.p=0.7,
#                                   pm.n=25, pm.p=0.3
#                      )) {
#   ## Extract control parameters:
#   default <- formals(sys.function())$control
#   control <- steady_state_emoa_control(f, lower, upper, ..., control=control, default=default)
#   control <- sbx_control(f, upper, lower, ..., control=control, default=default)
#   control <- pm_control(f, upper, lower, ..., control=control, default=default)  
#   control$ref <- emoa:::coalesce(control[["ref"]], rep(11, control$d))
#   
#   ## Tracking variables:
#   X <- matrix(0, nrow=control$n, ncol=control$maxeval)
#   Y <- matrix(0, nrow=control$d, ncol=control$maxeval)
#   dob <- rep(-1L, control$maxeval)
#   eol <- rep(-1L, control$maxeval)
#   
#   ## Random inital population:
#   X[, 1:control$mu] <- replicate(control$mu, runif(control$n, lower, upper))
#   Y[, 1:control$mu] <- sapply(1:control$mu, function(i) f(X[,i]))
#   
#   neval <- control$mu       ## Count the number of function evaluations
#   active <- 1:control$mu    ## Indices of individuals that are in the current pop.
#   
#   ## Save some common control parameters into the current
#   ## environment. This saves a few msec of execution time...
#   crossover <- control$crossover
#   mutate <- control$mutate
#   maxeval <- control$maxeval
#   logger <- control$logger
#   
#   logger$start("sms_emoa")
#   while(neval < maxeval) {
#     ############################################################
#     ## Variation:
#     parents <- sample(active, 2)
#     child <- crossover(X[, parents])[,sample(c(1, 2), 1)]
#     x <- mutate(child)
#     
#     ## Add new individual:
#     neval <- neval + 1
#     X[, neval] <- x
#     Y[, neval] <- f(x)
#     dob[neval] <- neval ## For a steady state emoa this is trivial...
#     active <- c(active, neval)
#     
#     ############################################################
#     ## Selection:
#     i <- nds_hv_selection(Y[, active])
#     
#     ## Remove the i-th active individual:
#     eol[active[i]] <- neval
#     active <- active[-i]
#     
#     ############################################################
#     ## Logging:    
#     logger$step()
#   }
#   logger$stop()
#   
#   res <- structure(list(X=X, Y=Y,
#                         dob=dob,
#                         eol=eol,
#                         par=X[,active], value=Y[,active]),
#                    class="emoa_result")
#   return(res)
# }
# 
# 
# zdt3 <- function (x) {
#   dim <- length(x)
#   y1 <- x[1]
#   g <- 1 + (9 * mean(x[2:dim]))
#   y2 <- g * (1 - sqrt(y1/g) - (y1/g) * sin(10 * pi * y1))
#   return(c(y1, y2))
# }
# 
# doit <- function() {
#   r1 <- sms_emoa(zdt3, rep(0, 5), rep(1, 5), control=control)
#   f1 <- nondominated_points(r1$value)
#   
#   
#   print(r1$value)
#   
#   
#   hv.sms <- dominated_hypervolume(f1, c(11, 11))
#   nopt.sms <- ncol(f1)
#   
#   r2 <- nsga2(zdt3, 5, 2,
#               lower.bounds=rep(0, 5), upper.bounds=rep(1, 5),
#               popsize=mu,
#               generations=gen)
#   f2 <- nondominated_points(t(r2$value))
#   
#   hv.nsga2 <- dominated_hypervolume(f2, c(11, 11))
#   nopt.nsga2 <- ncol(f2)
#   
#   data.frame(method=c("sms_emoa", "nsga2"),
#              hypervolume=c(hv.sms, hv.nsga2),
#              nopt=c(nopt.sms, nopt.nsga2))
# }
# 
# mu <- 100L
# gen <- 100L
# control <- list(mu=mu, maxeval=mu*gen,
#                 logger=emoa_console_logger(5000L))
# 
# res <- replicate(10L, doit(), simplify=FALSE)
# df <- Reduce(rbind, res)
# 
# plt <- qplot(hypervolume, data=subset(df, hypervolume > 127), fill=method, geom="density")
# plt




