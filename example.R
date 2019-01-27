library(plotly)
library(purrr)
source("./dtwCostMatrices.R")
source("./optimalWarpPath.R")

# Generate two time series
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

# Plot contour with (1, 1) in lower left
par(mfrow = c(1, 2))
image(x = 1:nrow(local.cost.m), 
      y = 1:ncol(local.cost.m), 
      z = local.cost.m, 
      col = gray.colors(15), 
      xlab = "i", ylab = "j", main = "Local cost matrix")
lines(map(path, 1), 
      map(path, 2), 
      type = "b", pch = 16, col = "black")
image(x = 1:nrow(acc.cost.m), 
      y = 1:ncol(acc.cost.m),
      z = acc.cost.m, 
      col = gray.colors(15), 
      xlab = "i", ylab = "j", main = "Accumulated cost matrix")
lines(map(path, 1), 
      map(path, 2), 
      type = "b", pch = 16, col = "black")
plot_ly(z = t(local.cost.m)) %>% add_surface()
plot_ly(z = t(acc.cost.m)) %>% add_surface()

# Plot stacked surface
plot_ly() %>% 
  add_surface(z = t(local.cost.m) - 1, 
              opacity = 1, 
              showscale = F) %>% 
  add_surface(x = 1:ncol(t(acc.cost.m)), 
              y = 1:nrow(t(acc.cost.m)), 
              z = t(acc.cost.m), 
              opacity = 1, 
              showscale = F) %>% 
  add_trace(x = map(path, 1), 
            y = map(path, 2), 
            z = (unlist(map(path, 3)) + 0.5), 
            type = "scatter3d", 
            mode = "lines", 
            line = list(color = "black"), 
            showlegend = F) %>% 
  layout(title = "Accumulated cost matrix")

# Plot contours
plot_ly() %>% 
  add_contour(x = 1:ncol(t(acc.cost.m)), 
              y = 1:nrow(t(acc.cost.m)), 
              z = t(acc.cost.m), 
              contours = list(showlabels = T), 
              showscale = F) %>% 
  add_trace(x = map(path, 1), 
            y = map(path, 2), 
            type = "scatter", mode = "lines+markers", 
            line = list(color = "black"), 
            marker = list(color = "black"), 
            showlegend = F) %>% 
  layout(title = "Accumulated cost matrix")
