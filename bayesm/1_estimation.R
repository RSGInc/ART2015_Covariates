###################################
# Setup and Data Preparation

# Set the working directory and load necessary packages
setwd("C://Users/jeff.keller/Desktop/ART2015_Covariates/bayesm/")
library(bayesm)
library(dummies)

# Fetch the choice data
load("../RSGHB/choicedata.RData")

# Remove hold-out sample
choicedata <- cd[out_hold == 0 & sequence != 3]
choicedata <- as.data.frame(choicedata)

# Create a choice variable
y <- as.numeric(choicedata$alt) * c(1, 0)[as.numeric(choicedata$choice)]
y <- aggregate(y, by = list(choicedata$sequence, choicedata$parent_id), sum)
y <- 1*(y[, 3] == 2)

# Create a matrix of independent variables
cbind(
     parent_id = choicedata[, "parent_id"],
     sequence = choicedata[, "sequence"],
     alt = choicedata[, "alt"],
     dummy(choicedata[, "pos"]),
     dummy(choicedata[, "att1"]),
     dummy(choicedata[, "att2"]),
     dummy(choicedata[, "att3"]),
     dummy(choicedata[, "att4"]),
     dummy(choicedata[, "att5"]),
     dummy(choicedata[, "att6"]),
     dummy(choicedata[, "att7"]),
     dummy(choicedata[, "att8"])
) -> design

# Organize data as binary logit
design[design[, "alt"] == 2, -c(1:3)] <- -1 * design[design[, "alt"] == 2, -c(1:3)]
design <- aggregate(design, by = list(design[, "sequence"], design[, "parent_id"]), sum)[, -c(1:2, 5)]
design[, 1:2] <- design[, 1:2] / 2 

# Create list of choice data for use with bayesm
lgtdata <- list()
counter <- 1
for(id in sort(unique(design[, "parent_id"]))) {
     lgtdata[[counter]] <- list(y = y[design[, "parent_id"] == id],
                                X = as.matrix(design[design[, "parent_id"] == id, -c(1:3, 5, 11, 16, 21, 25, 29, 33, 36)]))
     counter <- counter + 1
}

# Create matrix of covariates
Z <- matrix(0, 404, 2)
counter <- 1
for(id in sort(unique(design[, "parent_id"]))) {
  z <- as.numeric(choicedata$cov[1 + (counter - 1) * 7])
  Z[counter, ] <- 1 * c(z == 2, z == 3)
  counter <- counter + 1
}

###################################
### Model Controls and Settings
R <- 20000
keep <- 2

##############################
### Estimate Models

# Without Covariates
set.seed(1987)
NCV_bay <- rhierBinLogit(Data = list(lgtdata = lgtdata), Mcmc = list(R = R, keep = keep))

# With Covariates
set.seed(1987)
CV_bay <- rhierBinLogit(Data = list(lgtdata = lgtdata, Z = Z), Mcmc = list(R = R, keep = keep))

# Save Models
save(NCV_bay, file = "NCV_bay.RData")
save(CV_bay, file = "CV_bay.RData")

# # organize results
# cat("Summary of Delta draws", fill = TRUE)
# summary(model.nocovariate$Deltadraw)
# cat("Summary of Vbeta draws", fill = TRUE)
# summary(model.nocovariate$Vbetadraw)
# 
# model.nocovariate$beta
# 
# plot(model.nocovariate$Deltadraw)
# plot(model.nocovariate$betadraw)
# plot(model.nocovariate$Vbetadraw)
# 
# write.table(model.nocovariate$Deltadraw, "bayesm_nocovariate_delta.csv", sep = "," , row.names = FALSE)
# write.table(model.nocovariate$Vbetadraw, "bayesm_nocovariate_vbeta.csv", sep = "," , row.names = FALSE)
# 
# results <- matrix(0, 404, 28)
# 
# for(p in 1:404) {
#   results[p, ] <- c(sort(unique(choicedata$parent_id))[p], rowMeans(model.nocovariate$betadraw[p, 1:27, 5001:10000]))
# }
# write.table(results, "bayesm_nocovariate_beta.csv", sep = ",", row.names = FALSE)
# 
# write.table(model.covariate$Deltadraw, "bayesm_covariate_delta.csv", sep = ",", row.names = FALSE)
# write.table(model.covariate$Vbetadraw, "bayesm_covariate_vbeta.csv", sep = ",", row.names = FALSE)
# 
# results <- matrix(0, 404, 28)
# 
# for(p in 1:404) {
#   results[p, ] <- c(sort(unique(choicedata$parent_id))[p], rowMeans(model.covariate$betadraw[p, 1:27, 5001:10000]))
# }
# write.table(results, "bayesm_covariate_beta.csv", sep = ",", row.names = FALSE)