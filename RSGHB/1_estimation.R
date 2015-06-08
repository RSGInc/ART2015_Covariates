###################################
# Setup and Data Preparation

# Set the working directory and load necessary packages
setwd("C:/Users/jeff.keller/Desktop/ART2015_Covariates/RSGHB/")
library(RSGHB)
library(data.table)

# Fetch the choice data
load("choicedata.RData")

# Remove hold-out sample
choicedata <- cd[out_hold == 0 & sequence != 3]

###################################
# Defining the Likelihood Functions

source("0_likelihood_functions.R")

###################################
### Model Controls and Settings

# Random parameter names
params <- c("pos.left", "att1.Level1", "att1.Level2", "att1.Level3", "att1.Level4", 
            "att1.Level5", "att2.Level1", "att2.Level2", "att2.Level3", "att2.Level4",
            "att3.Level1", "att3.Level2", "att3.Level3", "att3.Level4", "att4.Level1",
            "att4.Level2", "att4.Level3", "att5.Level1", "att5.Level2", "att5.Level3",
            "att6.Level1", "att6.Level2", "att6.Level3", "att7.Level1", "att7.Level2",
            "att8.Level1", "att8.Level2")

# Covariate parameter names
covariates <- c("cov1_pos.left", "cov2_pos.left", "cov1_att1.Level1", "cov2_att1.Level1", "cov1_att1.Level2", "cov2_att1.Level2",
                "cov1_att1.Level3", "cov2_att1.Level3", "cov1_att1.Level4", "cov2_att1.Level4", "cov1_att1.Level5", "cov2_att1.Level5", "cov1_att2.Level1",
                "cov2_att2.Level1", "cov1_att2.Level2", "cov2_att2.Level2", "cov1_att2.Level3", "cov2_att2.Level3", "cov1_att2.Level4", "cov2_att2.Level4",
                "cov1_att3.Level1", "cov2_att3.Level1", "cov1_att3.Level2", "cov2_att3.Level2", "cov1_att3.Level3", "cov2_att3.Level3", "cov1_att3.Level4", "cov2_att3.Level4",
                "cov1_att4.Level1", "cov2_att4.Level1", "cov1_att4.Level2", "cov2_att4.Level2", "cov1_att4.Level3", "cov2_att4.Level3", "cov1_att5.Level1", "cov2_att5.Level1",
                "cov1_att5.Level2", "cov2_att5.Level2", "cov1_att5.Level3", "cov2_att5.Level3", "cov1_att6.Level1", "cov2_att6.Level1", "cov1_att6.Level2", "cov2_att6.Level2", "cov1_att6.Level3", "cov2_att6.Level3",
                "cov1_att7.Level1", "cov2_att7.Level1", "cov1_att7.Level2", "cov2_att7.Level2", "cov1_att8.Level1", "cov2_att8.Level1", "cov1_att8.Level2", "cov2_att8.Level2")

# ID variable represents the individuals associated with the output of the Likelihood functions
ID <- data.frame(unique(choicedata[, list(parent_id, sequence)]))
names(ID) <- c("ID", "sequence")

##############################
### Estimate Models

# Without Covariates (Default Settings) 4128 seconds
system.time(NCV_def <- doHB(Likelihoods_NCV, choicedata = ID,
                            control = list(gVarNamesFixed = NULL, gVarNamesNormal = params,
                                           modelname = "NCV_def", nodiagnostics = TRUE,
                                           gNCREP = 100000, gNEREP = 10000, gSeed = 1991)))

# Without Covariates (Modified Settings)
system.time(NCV_mod <- doHB(Likelihoods_NCV, choicedata = ID,
                            control = list(gVarNamesFixed = NULL, gVarNamesNormal = params,
                                           modelname = "NCV_mod", nodiagnostics = TRUE,
                                           gNCREP = 100000, gNEREP = 10000, gSeed = 1991,
                                           degreesOfFreedom = 0, priorVariance = 27)))

# With Covariates (Default Settings)
system.time(CV_def <- doHB(Likelihoods_CV, choicedata = ID,
                           control = list(gVarNamesFixed = covariates, gVarNamesNormal = params,
                                          modelname = "CV_def", nodiagnostics = TRUE,
                                          gNCREP = 100000, gNEREP = 10000, gSeed = 1991)))

# With Covariates (Modified Settings)
system.time(CV_mod <- doHB(Likelihoods_CV, choicedata = ID,
                           control = list(gVarNamesFixed = covariates, gVarNamesNormal = params,
                                          modelname = "CV_mod", nodiagnostics = TRUE,
                                          gNCREP = 100000, gNEREP = 10000, gSeed = 1991,
                                          degreesOfFreedom = 0, priorVariance = 27)))

# Save Models
save(NCV_def, file = "NCV_def.RData")
save(NCV_mod, file = "NCV_mod.RData")
save(CV_def,  file = "CV_def.RData")
save(CV_mod,  file = "CV_mod.RData")

rm(list = ls())