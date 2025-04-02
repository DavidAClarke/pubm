## Fit GLMs
fitglms <- function(zz, xx, xx_ppd, BB, nphi, nsill){ # xx = list of fine covariates
  
  be_est <- vector("list", length = 4)
  nn <- length(zz)
  KK <- ncol(xx_ppd[[1]][[1]][[1]])
  
  xxdf <- data.frame(zz = zz)
  for(d in 1:length(xx)){
    
    xxv <- xx[[d]][,nphi,nsill]
    xxdf <- data.frame(xxdf,xxv)
    
  }
  
  ns <- names(xxdf)
  y <- ns[1]
  x <- paste(ns[2:length(ns)], collapse = "+")
  f <- as.formula(paste(y,"~", x))
  
  #-------------------------------------------------------------------------------
  # GLM-oracle
  #-------------------------------------------------------------------------------
  # Fit a glm
  orig_out <- glm(f, data = xxdf, family = poisson)
  orig_be <- orig_out$coefficients
  
  # Parametric bootstrap
  orig_be_pb_sam <- array(NA, dim = c(length(orig_be), BB))
  xxs <- as.matrix(cbind(1, xxdf[,2:ncol(xxdf)]))
  orig_yy <- exp(xxs %*% orig_be)
  
  for (b in 1:BB) {
    zz_b <- rpois(nn, orig_yy)
    
    xxdf2 <- data.frame(zz_b, xxdf)
    ns2 <- names(xxdf2)
    y2 <- ns2[1]
    x2 <- paste(ns2[3:length(ns2)], collapse = "+")
    f2 <- as.formula(paste(y2,"~", x2))
    
    orig_out_b <- glm(f2, data = xxdf2, family = poisson)
    orig_be_pb_sam[, b] <- orig_out_b$coefficients
  }
  
  be_est[[1]] <- orig_be_pb_sam
  
  #-------------------------------------------------------------------------------
  # GLM-plugin
  #-------------------------------------------------------------------------------
  xx_plugdf <- data.frame(zz = zz)
  
  for(d in 1:length(xx_ppd)){
    
    xxv <- rowMeans(xx_ppd[[d]][[nphi]][[nsill]])
    xx_plugdf <- data.frame(xx_plugdf,xxv)
    
  }
  
  ns <- names(xx_plugdf)
  y <- ns[1]
  x <- paste(ns[2:length(ns)], collapse = "+")
  f <- as.formula(paste(y,"~", x))
  
  # Fit a glm
  plugin_out <- glm(f, data = xx_plugdf, family = poisson)
  plugin_be <- plugin_out$coefficients
  
  # Parametric bootstrap
  plugin_be_pb_sam <- array(NA, dim = c(length(plugin_be), BB))
  xx_plug <- as.matrix(cbind(1, xx_plugdf[,2:ncol(xx_plugdf)]))
  plugin_yy <- exp(xx_plug %*% plugin_be)
  
  for (b in 1:BB) {
    zz_b <- rpois(nn, plugin_yy)
    xx_plugdf2 <- data.frame(zz_b, xx_plugdf)
    ns <- names(xx_plugdf2)
    y2 <- ns[1]
    x2 <- paste(ns[3:length(ns)], collapse = "+")
    f2 <- as.formula(paste(y2,"~", x2))
    
    plugin_out_b <- glm(f2, data = xx_plugdf2, family = poisson)
    plugin_be_pb_sam[, b] <- plugin_out_b$coefficients
  }
  
  be_est[[2]] <- plugin_be_pb_sam
  
  #-------------------------------------------------------------------------------
  # GLM-ensemble
  #-------------------------------------------------------------------------------
  ppd_be_pb_sam_pool <- vector("list", length = KK)
  
  for (k in 1:KK) {
    
    xx_ensdf <- data.frame(zz = zz)
    
    for(d in 1:length(xx_ppd)){
      
      xxv <- xx_ppd[[d]][[nphi]][[nsill]][,k]
      xx_ensdf <- data.frame(xx_ensdf,xxv)
      
    }
    
    ns <- names(xx_ensdf)
    y <- ns[1]
    x <- paste(ns[2:length(ns)], collapse = "+")
    f <- as.formula(paste(y,"~", x))
    
    # Fit a glm
    ppd_out_k <- glm(f, data = xx_ensdf, family = poisson)
    ppd_be_k <- ppd_out_k$coefficients
    
    # Parametric bootstrap
    ppd_be_pb_sam_k <- array(NA, dim = c(length(ppd_be_k), BB))
    xx_ens <- as.matrix(cbind(1, xx_ensdf[,2:ncol(xx_ensdf)]))
    yy_k <- exp(xx_ens %*% ppd_be_k)
    
    for (b in 1:BB) {
      zz_b <- rpois(nn, yy_k)
      
      xx_ensdf2 <- data.frame(zz_b, xx_ensdf)
      ns <- names(xx_ensdf2)
      y2 <- ns[1]
      x2 <- paste(ns[3:length(ns)], collapse = "+")
      f2 <- as.formula(paste(y2,"~", x2))
      
      ppd_out_kb <- glm(f2, data = xx_ensdf2, family = poisson)
      ppd_be_pb_sam_k[, b] <- ppd_out_kb$coefficients
    }
    
    ppd_be_pb_sam_pool[[k]] <- ppd_be_pb_sam_k
    
  }
  
  be_est[[3]] <- ppd_be_pb_sam_pool
  
  #-------------------------------------------------------------------------------
  # GLM-Berkson
  #-------------------------------------------------------------------------------
  xx_bksdf <- data.frame(zz = zz)
  
  for(d in 1:length(xx_ppd)){
    
    xxv <- rowMeans(xx_ppd[[d]][[nphi]][[nsill]])
    xx_bksdf <- data.frame(xx_bksdf,xxv)
    
  }
  
  ns <- names(xx_bksdf)
  y <- ns[1]
  x <- paste(ns[2:length(ns)], collapse = "+")
  f <- as.formula(paste(y,"~", x))
  
  ppd_mu_be <- glm(f, data = xx_bksdf, family = poisson)$coefficients
  
  bks_be_pb_sam_pool <- vector("list", length = KK)
  
  for (k in 1:KK) {
    
    xx_bksdf2 <- data.frame(zz = zz)
    
    # Note: data.frame is subsequently zz, x1, delta_x1, x2, delta_x2, x3, delta_x3
    for(d in 1:length(xx_ppd)){
      
      xxv <- xx_ppd[[d]][[nphi]][[nsill]][,k]
      delta_k <- ppd_mu_be[d+1] * (xxv - xx_bksdf[d+1])
      xx_bksdf2 <- data.frame(xx_bksdf2,xxv,delta_k)
      
    }
    
    colnames(xx_bksdf2) <- c("zz", "x1", "delta_x1", "x2", "delta_x2", "x3", "delta_x3")
    ns <- names(xx_bksdf2)
    y <- ns[1]
    x <- paste(ns[c(2,4,6)], collapse = "+")
    f <- as.formula(paste(y,"~", x))
    
    uncert <- xx_bksdf2 %>% dplyr::select(delta_x1, delta_x2, delta_x3) %>% apply(1, sum)
    
    bks_out_k <- glm(f, data = xx_bksdf2, offset = uncert, family = poisson)
    bks_be_k <- bks_out_k$coefficients
    
    # Parametric bootstrap
    bks_be_pb_sam_k <- array(NA, dim = c(length(bks_be_k), BB))
    xx_bks <- as.matrix(cbind(1, xx_bksdf[,2:ncol(xx_bksdf)]))
    yy_k <- exp(xx_bks %*% bks_be_k + uncert)
    
    for (b in 1:BB) {
      zz_b <- rpois(nn, yy_k)
      
      xx_bksdf3 <- data.frame(zz_b, xx_bksdf2)
      ns <- names(xx_bksdf3)
      y <- ns[1]
      x <- paste(ns[c(3,5,7)], collapse = "+")
      f <- as.formula(paste(y,"~", x))
      
      bks_out_kb <- glm(f, data = xx_bksdf3, offset = uncert, family = poisson)
      bks_be_pb_sam_k[, b] <- bks_out_kb$coefficients
    }
    
    bks_be_pb_sam_pool[[k]] <- bks_be_pb_sam_k
    
  }
  
  be_est[[4]] <- bks_be_pb_sam_pool
  
  list(all_be = be_est, ppd_mu_be = ppd_mu_be)
  
}
