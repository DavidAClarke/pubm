#-------------------------------------------------------------------------------
# Load covariate data and response data
#-------------------------------------------------------------------------------

## Define file path
outpath <- here("data", "output")

## load true covariate
load(here(outpath, "covar_dat_full.rdata"))

## load downscaled covariate 
xx1_pred <- loadRData(here(outpath, "sim_all_x1_gpbau_pred.RData"))
xx2_pred <- loadRData(here(outpath, "sim_all_x2_gpbau_pred.RData"))
xx3_pred <- loadRData(here(outpath, "sim_all_x3_gpbau_pred.RData"))

## load response data
pois_rep <- loadRData(here(outpath, "pois_rep_data.rdata"))
#pois_rep <- loadRData(here(outpath, "pois_rep_data_all.rdata")) # multiple nphi and nsill
nrep <- length(pois_rep)

## combine covariate data into lists
#fine_list <- list(fine_xx_1, fine_xx_2, fine_xx_3)
xx_pred <- list(xx1_pred, xx2_pred, xx3_pred)

## Split training and testing data
fine_train_zz <- pois_rep[[1]]$train$fine_train_zz
fine_train_xx <- list(
  pois_rep[[1]]$train$fine_train_xx_1,
  pois_rep[[1]]$train$fine_train_xx_2,
  pois_rep[[1]]$train$fine_train_xx_3
)

#-------------------------------------------------------------------------------
# Fit GLMs
#-------------------------------------------------------------------------------
# Loop requires "pois_rep_data_all.rdata" to be loaded

pois_out_list_r <- list()

for(i in 1:2){
  
  pois_out_list_r[[i]] <- fitglms(zz = fine_train_zz, 
                         xx = fine_train_xx, 
                         xx_ppd = xx_pred, 
                         BB = 4, 
                         nphi = i, 
                         nsill = 1)
  
}



