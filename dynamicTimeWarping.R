# Dynamic time warping algorithm and calculation of optimal warp path

# Author: Michael Duprey
# Date: June 9, 2018

# Allows for comparison of two time series based on optimally time-shifted distance.

# Implementation based largely on Senin (2008):
# Senin, Pavel. "Dynamic time warping algorithm review."
# Information and Computer Science Department University of Hawaii at Manoa Honolulu, USA (2008).
# http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.465.4905&rep=rep1&type=pdf


library(plotly)
library(purrr)


# Dynamic time warping ----------------------------------------------------

dtwCostMatrices <- function(s, t) {
  m <- length(s)
  n <- length(t)
  
  # Define a distance/cost function
  d <- function(x, y) abs(x - y) # Euclidean distance
  # d <- function(x, y) (x - y)^2 # Squared Euclidean distance
  
  # Accumulated cost matrix
  acc.cost.m <- matrix(nrow = m, ncol = n)
  
  # Local cost matrix
  local.cost.m <- matrix(nrow = m, ncol = n)
  
  # Set initial boundary distance
  acc.cost.m[2:m, 1] <- Inf
  acc.cost.m[1, 2:n] <- Inf
  acc.cost.m[1, 1] <- 0
  
  # Alternative method of setting boundary distance
  # for(i in 2:m) {
  #   acc.cost.m[i, 1] <- acc.cost.m[i - 1, 1] + d(i, 1)
  # }
  # 
  # for(j in 2:n) {
  #   acc.cost.m[1, j] <- acc.cost.m[1, j - 1] + d(1, j)
  # }
  
  for(i in 2:m) {
    for(j in 2:n) {
      cost <- d(s[i], t[j])
      
      acc.cost.m[i, j] <- cost + min(acc.cost.m[i - 1, j], acc.cost.m[i, j - 1], acc.cost.m[i - 1, j - 1])
      local.cost.m[i, j] <- cost
    }
  }
  
  # Trim initial boundary distance and extract the minimal distance between s- and t-sequences
  dtw <- list("acc.cost.m" = acc.cost.m[2:m, 2:n], 
              "local.cost.m" = local.cost.m[2:m, 2:n], 
              "min.dist" = acc.cost.m[m, n])
  
  return(dtw)
}

optimalWarpPath <- function(acc.cost.m) {
  i <- nrow(acc.cost.m)
  j <- ncol(acc.cost.m)
  
  # Optimal warping path
  path <- list(c(i, j, acc.cost.m[i, j]))

  while(i > 1 & j > 1) {
    if(i == 2) {
      j <- j - 1
    } else if(j == 2) {
      i <- i - 1
    } else {
      if(acc.cost.m[i - 1, j] == min(acc.cost.m[i - 1, j], acc.cost.m[i, j - 1], acc.cost.m[i - 1, j - 1])) {
        i <- i - 1
      } else if(acc.cost.m[i, j - 1] == min(acc.cost.m[i - 1, j], acc.cost.m[i, j - 1], acc.cost.m[i - 1, j - 1])) {
        j <- j - 1
      } else {
        i <- i - 1
        j <- j - 1
      }
    }
    path <- c(path, list(c(i, j, acc.cost.m[i, j])))
  }
  path <- c(path, list(c(1, 1, 1)))
  
  return(path)
}

# Example usage
s <- sin(1:12)
t <- cos(1:12)

plot(s, type = "b", xlab = "Time", ylab = "Value", main = "s-sequence", xlim = c(0, 12), ylim = c(-5, 5))
plot(t, type = "b", xlab = "Time", ylab = "Value", main = "t-sequence", xlim = c(0, 12), ylim = c(-5, 5))

# Apply DTW over s- and t-sequences
dtw <- dtwCostMatrices(s, t)

# Extract cost matrices
acc.cost.m <- dtw$acc.cost.m
local.cost.m <- dtw$local.cost.m

# Calculate optimal warp path
path <- optimalWarpPath(acc.cost.m)


# Graphing ----------------------------------------------------------------

# Plot contour with (1, 1) in lower left
par(mfrow = c(1, 2))
image(x = 1:nrow(local.cost.m), y = 1:ncol(local.cost.m), z = local.cost.m, col = gray.colors(15), 
      xlab = "i", ylab = "j", main = "Local cost matrix")
lines(map(path, 1), map(path, 2), type = "b", pch = 16, col = "black")
image(x = 1:nrow(acc.cost.m), y = 1:ncol(acc.cost.m), z = acc.cost.m, col = gray.colors(15), 
      xlab = "i", ylab = "j", main = "Accumulated cost matrix")
lines(map(path, 1), map(path, 2), type = "b", pch = 16, col = "black")
plot_ly(z = t(local.cost.m)) %>% add_surface()
plot_ly(z = t(acc.cost.m)) %>% add_surface()

# Plot stacked surface
plot_ly() %>% 
  add_surface(z = t(local.cost.m) - 1, opacity = 1, showscale = F) %>% 
  add_surface(x = 1:ncol(t(acc.cost.m)), y = 1:nrow(t(acc.cost.m)), z = t(acc.cost.m), opacity = 1, showscale = F) %>% 
  add_trace(x = map(path, 1), y = map(path, 2), z = (unlist(map(path, 3)) + 0.5), 
            type = "scatter3d", mode = "lines", line = list(color = "black"), showlegend = F) %>% 
  layout(title = "Accumulated cost matrix")

# Plot contours
plot_ly() %>% 
  add_contour(x = 1:ncol(t(acc.cost.m)), y = 1:nrow(t(acc.cost.m)), z = t(acc.cost.m), contours = list(showlabels = T), showscale = F) %>% 
  add_trace(x = map(path, 1), y = map(path, 2), type = "scatter", mode = "lines+markers", 
            line = list(color = "black"), marker = list(color = "black"), showlegend = F) %>% 
  layout(title = "Accumulated cost matrix")
