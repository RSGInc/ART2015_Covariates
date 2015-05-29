# Without Covariates
Likelihoods_NCV <- function(fc, b) {  
  
  # Define Parameters
  cc <- 1
  beta_pos.left <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att1.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att1.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att1.Level3 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att1.Level4 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att1.Level5 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att2.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att2.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att2.Level3 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att2.Level4 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att3.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att3.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att3.Level3 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att3.Level4 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att4.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att4.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att4.Level3 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att5.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att5.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att5.Level3 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att6.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att6.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att6.Level3 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att7.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att7.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att8.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att8.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  
  # Calculate Utility
  choicedata[, U :=  beta_pos.left * (pos == "left") +
               beta_att1.Level1 * (att1 == "Level1") +
               beta_att1.Level2 * (att1 == "Level2") +
               beta_att1.Level3 * (att1 == "Level3") +
               beta_att1.Level4 * (att1 == "Level4") +
               beta_att1.Level5 * (att1 == "Level5") +
               beta_att2.Level1 * (att2 == "Level1") +
               beta_att2.Level2 * (att2 == "Level2") +
               beta_att2.Level3 * (att2 == "Level3") +
               beta_att2.Level4 * (att2 == "Level4") +
               beta_att3.Level1 * (att3 == "Level1") +
               beta_att3.Level2 * (att3 == "Level2") +
               beta_att3.Level3 * (att3 == "Level3") +
               beta_att3.Level4 * (att3 == "Level4") +
               beta_att4.Level1 * (att4 == "Level1") +
               beta_att4.Level2 * (att4 == "Level2") +
               beta_att4.Level3 * (att4 == "Level3") +
               beta_att5.Level1 * (att5 == "Level1") +
               beta_att5.Level2 * (att5 == "Level2") +
               beta_att5.Level3 * (att5 == "Level3") +
               beta_att6.Level1 * (att6 == "Level1") +
               beta_att6.Level2 * (att6 == "Level2") +
               beta_att6.Level3 * (att6 == "Level3") +
               beta_att7.Level1 * (att7 == "Level1") +
               beta_att7.Level2 * (att7 == "Level2") +
               beta_att8.Level1 * (att8 == "Level1") +
               beta_att8.Level2 * (att8 == "Level2")]
  
  # Calculate Probabilities
  choicedata[, P := exp(U)/sum(exp(U)), by = list(parent_id, sequence)]
  return(choicedata[choice == 1, P])
  
}

# With Covariates
Likelihoods_CV <- function(fc, b) {  
  
  # Define Fixed Parameters
  cc <- 1
  cov1_pos.left <- fc[cc]; cc <- cc + 1
  cov2_pos.left <- fc[cc]; cc <- cc + 1
  cov1_att1.Level1 <- fc[cc]; cc <- cc + 1
  cov2_att1.Level1 <- fc[cc]; cc <- cc + 1
  cov1_att1.Level2 <- fc[cc]; cc <- cc + 1
  cov2_att1.Level2 <- fc[cc]; cc <- cc + 1
  cov1_att1.Level3 <- fc[cc]; cc <- cc + 1
  cov2_att1.Level3 <- fc[cc]; cc <- cc + 1
  cov1_att1.Level4 <- fc[cc]; cc <- cc + 1
  cov2_att1.Level4 <- fc[cc]; cc <- cc + 1
  cov1_att1.Level5 <- fc[cc]; cc <- cc + 1
  cov2_att1.Level5 <- fc[cc]; cc <- cc + 1
  cov1_att2.Level1 <- fc[cc]; cc <- cc + 1
  cov2_att2.Level1 <- fc[cc]; cc <- cc + 1
  cov1_att2.Level2 <- fc[cc]; cc <- cc + 1
  cov2_att2.Level2 <- fc[cc]; cc <- cc + 1
  cov1_att2.Level3 <- fc[cc]; cc <- cc + 1
  cov2_att2.Level3 <- fc[cc]; cc <- cc + 1 
  cov1_att2.Level4 <- fc[cc]; cc <- cc + 1
  cov2_att2.Level4 <- fc[cc]; cc <- cc + 1
  cov1_att3.Level1 <- fc[cc]; cc <- cc + 1
  cov2_att3.Level1 <- fc[cc]; cc <- cc + 1
  cov1_att3.Level2 <- fc[cc]; cc <- cc + 1
  cov2_att3.Level2 <- fc[cc]; cc <- cc + 1
  cov1_att3.Level3 <- fc[cc]; cc <- cc + 1
  cov2_att3.Level3 <- fc[cc]; cc <- cc + 1
  cov1_att3.Level4 <- fc[cc]; cc <- cc + 1
  cov2_att3.Level4 <- fc[cc]; cc <- cc + 1
  cov1_att4.Level1 <- fc[cc]; cc <- cc + 1
  cov2_att4.Level1 <- fc[cc]; cc <- cc + 1
  cov1_att4.Level2 <- fc[cc]; cc <- cc + 1
  cov2_att4.Level2 <- fc[cc]; cc <- cc + 1
  cov1_att4.Level3 <- fc[cc]; cc <- cc + 1
  cov2_att4.Level3 <- fc[cc]; cc <- cc + 1
  cov1_att5.Level1 <- fc[cc]; cc <- cc + 1
  cov2_att5.Level1 <- fc[cc]; cc <- cc + 1
  cov1_att5.Level2 <- fc[cc]; cc <- cc + 1
  cov2_att5.Level2 <- fc[cc]; cc <- cc + 1
  cov1_att5.Level3 <- fc[cc]; cc <- cc + 1
  cov2_att5.Level3 <- fc[cc]; cc <- cc + 1 
  cov1_att6.Level1 <- fc[cc]; cc <- cc + 1
  cov2_att6.Level1 <- fc[cc]; cc <- cc + 1
  cov1_att6.Level2 <- fc[cc]; cc <- cc + 1
  cov2_att6.Level2 <- fc[cc]; cc <- cc + 1
  cov1_att6.Level3 <- fc[cc]; cc <- cc + 1
  cov2_att6.Level3 <- fc[cc]; cc <- cc + 1
  cov1_att7.Level1 <- fc[cc]; cc <- cc + 1
  cov2_att7.Level1 <- fc[cc]; cc <- cc + 1
  cov1_att7.Level2 <- fc[cc]; cc <- cc + 1
  cov2_att7.Level2 <- fc[cc]; cc <- cc + 1
  cov1_att8.Level1 <- fc[cc]; cc <- cc + 1
  cov2_att8.Level1 <- fc[cc]; cc <- cc + 1
  cov1_att8.Level2 <- fc[cc]; cc <- cc + 1
  cov2_att8.Level2 <- fc[cc]; cc <- cc + 1
  
  # Define Random Parameters
  cc <- 1
  beta_pos.left    <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att1.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att1.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att1.Level3 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att1.Level4 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att1.Level5 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att2.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att2.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att2.Level3 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att2.Level4 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att3.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att3.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att3.Level3 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att3.Level4 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att4.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att4.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att4.Level3 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att5.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att5.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att5.Level3 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att6.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att6.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att6.Level3 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att7.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att7.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att8.Level1 <- rep(b[,cc], each = 2); cc <- cc + 1
  beta_att8.Level2 <- rep(b[,cc], each = 2); cc <- cc + 1
  
  # Calculate Utility
  choicedata[, U :=  (beta_pos.left + cov1_pos.left * (cov == "cov1") + cov2_pos.left * (cov == "cov2")) * (pos == "left") +
               (beta_att1.Level1 + cov1_att1.Level1 * (cov == "cov1") + cov2_att1.Level1 * (cov == "cov2")) * (att1 == "Level1") +
               (beta_att1.Level2 + cov1_att1.Level2 * (cov == "cov1") + cov2_att1.Level2 * (cov == "cov2")) * (att1 == "Level2") +
               (beta_att1.Level3 + cov1_att1.Level3 * (cov == "cov1") + cov2_att1.Level3 * (cov == "cov2")) * (att1 == "Level3") +
               (beta_att1.Level4 + cov1_att1.Level4 * (cov == "cov1") + cov2_att1.Level4 * (cov == "cov2")) * (att1 == "Level4") +
               (beta_att1.Level5 + cov1_att1.Level5 * (cov == "cov1") + cov2_att1.Level5 * (cov == "cov2")) * (att1 == "Level5") +
               (beta_att2.Level1 + cov1_att2.Level1 * (cov == "cov1") + cov2_att2.Level1 * (cov == "cov2")) * (att2 == "Level1") +
               (beta_att2.Level2 + cov1_att2.Level2 * (cov == "cov1") + cov2_att2.Level2 * (cov == "cov2")) * (att2 == "Level2") +
               (beta_att2.Level3 + cov1_att2.Level3 * (cov == "cov1") + cov2_att2.Level3 * (cov == "cov2")) * (att2 == "Level3") +
               (beta_att2.Level4 + cov1_att2.Level4 * (cov == "cov1") + cov2_att2.Level4 * (cov == "cov2")) * (att2 == "Level4") +
               (beta_att3.Level1 + cov1_att3.Level1 * (cov == "cov1") + cov2_att3.Level1 * (cov == "cov2")) * (att3 == "Level1") +
               (beta_att3.Level2 + cov1_att3.Level2 * (cov == "cov1") + cov2_att3.Level2 * (cov == "cov2")) * (att3 == "Level2") +
               (beta_att3.Level3 + cov1_att3.Level3 * (cov == "cov1") + cov2_att3.Level3 * (cov == "cov2")) * (att3 == "Level3") +
               (beta_att3.Level4 + cov1_att3.Level4 * (cov == "cov1") + cov2_att3.Level4 * (cov == "cov2")) * (att3 == "Level4") +
               (beta_att4.Level1 + cov1_att4.Level1 * (cov == "cov1") + cov2_att4.Level1 * (cov == "cov2")) * (att4 == "Level1") +
               (beta_att4.Level2 + cov1_att4.Level2 * (cov == "cov1") + cov2_att4.Level2 * (cov == "cov2")) * (att4 == "Level2") +
               (beta_att4.Level3 + cov1_att4.Level3 * (cov == "cov1") + cov2_att4.Level3 * (cov == "cov2")) * (att4 == "Level3") +
               (beta_att5.Level1 + cov1_att5.Level1 * (cov == "cov1") + cov2_att5.Level1 * (cov == "cov2")) * (att5 == "Level1") +
               (beta_att5.Level2 + cov1_att5.Level2 * (cov == "cov1") + cov2_att5.Level2 * (cov == "cov2")) * (att5 == "Level2") +
               (beta_att5.Level3 + cov1_att5.Level3 * (cov == "cov1") + cov2_att5.Level3 * (cov == "cov2")) * (att5 == "Level3") +
               (beta_att6.Level1 + cov1_att6.Level1 * (cov == "cov1") + cov2_att6.Level1 * (cov == "cov2")) * (att6 == "Level1") +
               (beta_att6.Level2 + cov1_att6.Level2 * (cov == "cov1") + cov2_att6.Level2 * (cov == "cov2")) * (att6 == "Level2") +
               (beta_att6.Level3 + cov1_att6.Level3 * (cov == "cov1") + cov2_att6.Level3 * (cov == "cov2")) * (att6 == "Level3") +
               (beta_att7.Level1 + cov1_att7.Level1 * (cov == "cov1") + cov2_att7.Level1 * (cov == "cov2")) * (att7 == "Level1") +
               (beta_att7.Level2 + cov1_att7.Level2 * (cov == "cov1") + cov2_att7.Level2 * (cov == "cov2")) * (att7 == "Level2") +
               (beta_att8.Level1 + cov1_att8.Level1 * (cov == "cov1") + cov2_att8.Level1 * (cov == "cov2")) * (att8 == "Level1") +
               (beta_att8.Level2 + cov1_att8.Level2 * (cov == "cov1") + cov2_att8.Level2 * (cov == "cov2")) * (att8 == "Level2")]
  
  # Calculate Probabilities
  choicedata[, P := exp(U)/sum(exp(U)), by = list(parent_id, sequence)]
  return(choicedata[choice == 1, P])
  
}
