## Simulate fine-resolution data
sim_vals <- function(fine_grid, tran = F, fun, range_vals, sill_vals){
  
  nn <- nrow(fine_grid)
  
  fine_xx_1 <- array(NA, dim = c(nn, length(range_vals), length(sill_vals)))
  
  for (i in seq(range_vals)) {
    
    for (j in seq(sill_vals)) {
      
      if(tran == T){
        
        xx <- t(SpatialExtremes::rgp(1, fine_grid, cov.mod = "powexp", nugget = 0, smooth = 1,
                                     sill = sill_vals[j], range = range_vals[i]))
        fine_xx_1[, i, j] <- fun(xx)
        
      }
      if(tran == F){
        
        xx <- t(SpatialExtremes::rgp(1, fine_grid, cov.mod = "powexp", nugget = 0, smooth = 1,
                                     sill = sill_vals[j], range = range_vals[i]))
        fine_xx_1[, i, j] <- xx
        
      }
    }
  }
  return(fine_xx_1)
}


## Aggregate simulated fine-resolution data
coar_vals <- function(fine_vals, tran = F, fun, range_vals, sill_vals){
  
  coar_xx <- array(NA, dim = c(mm, length(range_vals), length(sill_vals)))
  
  for (i in seq(range_vals)) {
    
    for (j in seq(sill_vals)) {
      
      if(tran == T){
        
        xx <- fine_vals[, i, j]
        
        xxb <- aggregate(xx, by = list(fine_grid$par_label), mean)$x
        
        coar_xx[, i, j] <- fun(xxb)
        
      }
      if(tran == F){
        
        xx <- fine_vals[, i, j]
        
        xxb <- aggregate(xx, by = list(fine_grid$par_label), mean)$x
        
        coar_xx[, i, j] <- xxb
        
      }
    }
  }
  return(coar_xx)
}
