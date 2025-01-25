

plot_predicteds<-function(oxymat, pred_ypred){
  
  #turn preds into a matrix
  indices <- gsub("ypred\\[|\\]", "", names(pred_ypred))
  index_matrix <- do.call(rbind, strsplit(indices, ","))
  
  # Create a matrix
  row_indices <- as.numeric(index_matrix[, 1])
  col_indices <- as.numeric(index_matrix[, 2])
  
  oxypred_matrix <- matrix(NA, nrow = max(row_indices), ncol = max(col_indices))
  for (i in seq_along(pred_ypred)) {
    oxypred_matrix[row_indices[i], col_indices[i]] <- pred_ypred[i]
  }
  
  
  
  par(mfrow=c(3,10), mai=c(0.2,0.2,0.1,0.1), mgp=c(2,1,0))   
  for (d in 1:nday){
    
    plot(seq(1:144), oxymat[,d] )
         lines(seq(1:144), oxypred_matrix[,d] )
    
  }
  
}


plot_predicteds(oxymat=oxymat, pred_ypred=pred_ypred)



plot_resids<-function(oxymat, pred_mu){
  
  #turn preds into a matrix
  indices <- gsub("mu\\[|\\]", "", names(pred_mu))
  index_matrix <- do.call(rbind, strsplit(indices, ","))
  
  # Create a matrix
  row_indices <- as.numeric(index_matrix[, 1])
  col_indices <- as.numeric(index_matrix[, 2])
  
  oxypred_matrix <- matrix(NA, nrow = max(row_indices), ncol = max(col_indices))
  for (i in seq_along(pred_mu)) {
    oxypred_matrix[row_indices[i], col_indices[i]] <- pred_mu[i]
  }
  
  
  
  par(mfrow=c(3,10), mai=c(0.2,0.2,0.1,0.1), mgp=c(2,1,0))   
  for (d in 1:nday){
    eta<- oxymat[,d]-oxypred_matrix[,d]
    plot(seq(1:144), eta, type="o" )
    
    
  }
  
}

plot_resids(oxymat=oxymat, pred_mu=mu_summary[,1])



plot_resids_withlight<-function(oxymat, pred_mu, light_mat){
  
  #turn preds into a matrix
  indices <- gsub("mu\\[|\\]", "", names(pred_mu))
  index_matrix <- do.call(rbind, strsplit(indices, ","))
  
  # Create a matrix
  row_indices <- as.numeric(index_matrix[, 1])
  col_indices <- as.numeric(index_matrix[, 2])
  
  oxypred_matrix <- matrix(NA, nrow = max(row_indices), ncol = max(col_indices))
  for (i in seq_along(pred_mu)) {
    oxypred_matrix[row_indices[i], col_indices[i]] <- pred_mu[i]
  }
  
  
  
  par(mfrow=c(3,10), mai=c(0.2,0.2,0.1,0.1), mgp=c(2,1,0))   
  for (d in 1:nday){
    eta<- oxymat[,d]-oxypred_matrix[,d]
    plot(light_mat[,d], abs(eta) )
    
    
  }
  
}

plot_resids_withlight(oxymat=oxymat, pred_mu=mu_summary[,1], light_mat = matrix(oxy_s$light, nrow=ntime) )


##Function to simply get matix of mu
mu_matrix<- function(pred_mu){
  
  #turn preds into a matrix
  indices <- gsub("mu\\[|\\]", "", names(pred_mu))
  index_matrix <- do.call(rbind, strsplit(indices, ","))
  
  # Create a matrix
  row_indices <- as.numeric(index_matrix[, 1])
  col_indices <- as.numeric(index_matrix[, 2])
  
  oxypred_matrix <- matrix(NA, nrow = max(row_indices), ncol = max(col_indices))
  for (i in seq_along(pred_mu)) {
    oxypred_matrix[row_indices[i], col_indices[i]] <- pred_mu[i]
  }
  oxypred_matrix
  
}

mu_predict_matrix<-mu_matrix(pred_mu=pred_mu)

mu_summary <- summary(oipi_fit, pars = c("mu"), probs = c(0.5))$summary




resids<-oxymat-mu_predict_matrix

plot(matrix(oxy_s$light, nrow=ntime), log(abs(resids)) )


plot(seq(1:144),oxymat[,1])
lines(seq(1:144), oxypred_matrix [,1])

plot(seq(1:143), diff(oxymat[,1]))
lines(seq(1:143), diff(oxypred_matrix [,1]))

plot(seq(1:144),oxymat[,1]-oxypred_matrix [,1])

 