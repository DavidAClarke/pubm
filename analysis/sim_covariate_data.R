#-------------------------------------------------------------------------------
# Simulate covariate data
#-------------------------------------------------------------------------------
pkgs <- c("SpatialExtremes", "ggplot2", "dplyr", "here")
lapply(pkgs, require, character.only = T)

outpath <- here("data", "output")

set.seed(42)

#-------------------------------------------------------------------------------
# Generate true fine-resolution covariate fields
#-------------------------------------------------------------------------------
range_vals <- seq(0.1, 0.8, by = 0.1)
sill_vals <- seq(1, 5, by = 1)

x_seq <- seq(0, 1, length = 50)
y_seq <- seq(0, 1, length = 60)

fine_grid <- as.matrix(expand.grid(x_seq, y_seq))
colnames(fine_grid) <- c("lon", "lat")
nn <- nrow(fine_grid)

source("R/sim_covar_data.R")

fine_xx_1 <- sim_vals(fine_grid = fine_grid, 
                      tran = F, 
                      range_vals = range_vals, 
                      sill_vals = sill_vals)

fine_xx_2 <- sim_vals(fine_grid = fine_grid, 
                      tran = T, 
                      fun = exp, # values > 0
                      range_vals = range_vals, 
                      sill_vals = sill_vals)

fine_xx_3 <- sim_vals(fine_grid = fine_grid, 
                      tran = T, 
                      fun = plogis, # values [0,1]
                      range_vals = range_vals, 
                      sill_vals = sill_vals)

#-------------------------------------------------------------------------------
# Aggregate to coarse-resolution data
#-------------------------------------------------------------------------------
kk_x <- 5
kk_y <- 6
mm <- kk_x * kk_y

fine_grid <- as.data.frame(fine_grid)
fine_grid$x_cut <- as.numeric(cut(fine_grid[, "lon"], kk_x))
fine_grid$y_cut <- as.numeric(cut(fine_grid[, "lat"], kk_y))
grid_index <- expand.grid(x_cut = 1:kk_x, y_cut = 1:kk_y)
grid_index$par_label <- as.numeric(row.names(grid_index))
fine_grid <- left_join(fine_grid, grid_index)

coar_lon <- aggregate(fine_grid$lon, by = list(fine_grid$par_label), mean)$x
coar_lat <- aggregate(fine_grid$lat, by = list(fine_grid$par_label), mean)$x
coar_grid <- data.frame(lon = coar_lon, lat = coar_lat)

coar_xx_1 <- coar_vals(fine_vals = fine_xx_1, 
                       tran = F, 
                       range_vals = range_vals, 
                       sill_vals = sill_vals)

coar_xx_2 <- coar_vals(fine_vals = fine_xx_2, 
                       tran = T, 
                       fun = log, # back transform from > 0
                       range_vals = range_vals, 
                       sill_vals = sill_vals)

coar_xx_3 <- coar_vals(fine_vals = fine_xx_3, 
                       tran = T, 
                       fun = qlogis, # back transform from [0,1]
                       range_vals = range_vals, 
                       sill_vals = sill_vals)

#-------------------------------------------------------------------------------
# Output covariate data
#-------------------------------------------------------------------------------
save(fine_grid, coar_grid, 
     fine_xx_1,fine_xx_2,fine_xx_3,
     coar_xx_1,coar_xx_2,coar_xx_3,
     range_vals, sill_vals,
     file = here(outpath, "covar_dat_full.rdata"))