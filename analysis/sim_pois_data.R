#-------------------------------------------------------------------------------
# Simulating (Poisson) response data
#-------------------------------------------------------------------------------
inpath <- here("data", "output")
outpath <- here("data", "output")

set.seed(42)

#-------------------------------------------------------------------------------
# Load covariate data
#-------------------------------------------------------------------------------
load(here(inpath, "covar_dat_full.rdata"))

nphi <- dim(fine_xx_1)[2]
nsill <- dim(fine_xx_1)[3]
#-------------------------------------------------------------------------------
# Generate poisson data
# Split the data into training and testing set
#-------------------------------------------------------------------------------
b0 <- 1
b1 <- 0.5
b2 <- 0.5
b3 <- 0.5
true_be <- c(b0, b1, b2, b3) # Values don't matter at this point.
nn <- nrow(fine_xx_1)
nrep <- 100
rep_data <- vector("list", length = nrep)
train_idx <- sample(1:nn, nn / 5)

for (r in 1:nrep) {
  
  fine_yy <- array(NA, dim = c(nn, nphi, nsill))
  fine_zz <- array(NA, dim = c(nn, nphi, nsill))
  
  for (i in 1:nphi) {
    
    for (j in 1:nsill) {
      
      xx1 <- fine_xx_1[, i, j]
      xx2 <- fine_xx_2[, i, j]
      xx3 <- fine_xx_3[, i, j]
      yy <- exp(cbind(1, xx1,xx2,xx3) %*% true_be)
      zz <- rpois(nn, yy)
      
      fine_yy[, i, j] <- yy
      fine_zz[, i, j] <- zz
      
    }
    
  }
  
  all <- list(fine_xx_1 = fine_xx_1, 
              fine_xx_2 = fine_xx_2,
              fine_xx_3 = fine_xx_3,
              fine_yy = fine_yy, 
              fine_zz = fine_zz)
  
  train <- list(fine_train_xx_1 = fine_xx_1[train_idx, , ],
                fine_train_xx_2 = fine_xx_2[train_idx, , ],
                fine_train_xx_3 = fine_xx_3[train_idx, , ],
                fine_train_yy = fine_yy[train_idx, , ], 
                fine_train_zz = fine_zz[train_idx, , ])
  
  test <- list(fine_test_xx_1 = fine_xx_1[-train_idx, , ], 
               fine_test_xx_2 = fine_xx_2[-train_idx, , ], 
               fine_test_xx_3 = fine_xx_3[-train_idx, , ], 
               fine_test_yy = fine_yy[-train_idx, , ],
               fine_test_zz = fine_zz[-train_idx, , ])
  
  rep_data[[r]] <- list(all = all, train = train, test = test,
                        train_idx = train_idx, 
                        true_be = true_be)
  
}

save(rep_data, file = here(outpath, "pois_rep_data.rdata"))
