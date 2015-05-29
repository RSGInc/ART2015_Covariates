##############################################
# Load packages, data, and model objects

setwd("C:/Users/jeff.keller/Desktop/ART2015_Covariates/RSGHB/")

library(RSGHB)
library(data.table)
library(MASS)

load("CV_def.RData")
load("CV_mod.RData")
load("NCV_def.RData")
load("NCV_mod.RData")
load("choicedata.RData")

source("0_likelihood_functions.R")

##############################################
# Model Summaries

# Lower-level model RLH
RLH_table <- data.frame(RLH = c(mean(CV_def [["C"]][, "RLH"]),
                                mean(NCV_def[["C"]][, "RLH"]),
                                mean(CV_mod [["C"]][, "RLH"]),
                                mean(NCV_mod[["C"]][, "RLH"])),
                        row.names = c("CV_def", "NCV_def", "CV_mod", "NCV_mod"))

##############################################
# Upper-level model simulated LL

N <- 404
nDraws <- 1000
choicedata <- cd[out_hold == 0 & sequence != 3]
LLSimTable <- data.frame(AvgLL = rep(NA, 4), row.names = c("CV_def", "NCV_def", "CV_mod", "NCV_mod"))

# With Covariates (default settings)
ll <- rep(NA, nDraws)
fc <- colMeans(CV_def[["F"]][, -1])
A <- colMeans(CV_def[["A"]][, -1])
D <- apply(X = CV_def[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
  b <- mvrnorm(N, mu = A, Sigma = D)
  b <- b[rep(1:nrow(b), each = 7), ]
  ll[draw] <- sum(log(Likelihoods_CV(fc = fc, b = b)))
} # This might take a few minutes

LLSimTable["CV_def", "AvgLL"] <- mean(ll)

# Without Covariates (default settings)
ll <- rep(NA, nDraws)
fc <- NULL
A <- colMeans(NCV_def[["A"]][, -1])
D <- apply(X = NCV_def[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
  b <- mvrnorm(N, mu = A, Sigma = D)
  b <- b[rep(1:nrow(b), each = 7), ]
  ll[draw] <- sum(log(Likelihoods_NCV(fc = fc, b = b)))
} # This might take a few minutes

LLSimTable["NCV_def", "AvgLL"] <- mean(ll)

# With Covariates (modified settings)
ll <- rep(NA, nDraws)
fc <- colMeans(CV_mod[["F"]][, -1])
A <- colMeans(CV_mod[["A"]][, -1])
D <- apply(X = CV_mod[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
  b <- mvrnorm(N, mu = A, Sigma = D)
  b <- b[rep(1:nrow(b), each = 7), ]
  ll[draw] <- sum(log(Likelihoods_CV(fc = fc, b = b)))
} # This might take a few minutes

LLSimTable["CV_mod", "AvgLL"] <- mean(ll)

# Without Covariates (modified settings)
ll <- rep(NA, nDraws)
fc <- NULL
A <- colMeans(NCV_mod[["A"]][, -1])
D <- apply(X = NCV_mod[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
  b <- mvrnorm(N, mu = A, Sigma = D)
  b <- b[rep(1:nrow(b), each = 7), ]
  ll[draw] <- sum(log(Likelihoods_NCV(fc = fc, b = b)))
} # This might take a few minutes

LLSimTable["NCV_mod", "AvgLL"] <- mean(ll)

##############################################
# Complete Respondent Hold-Out Sample

N <- 45
nDraws <- 1000
choicedata <- cd[out_hold == 1]
RespObsTable <- data.frame(HitRate = rep(NA, 4),
                           AvgP = rep(NA, 4),
                           row.names = c("CV_def", "NCV_def", "CV_mod", "NCV_mod"))

# With Covariates (default settings)
hitrate <- rep(NA, nDraws)
prob <- rep(NA, nDraws)
fc <- colMeans(CV_def[["F"]][, -1])
A <- colMeans(CV_def[["A"]][, -1])
D <- apply(X = CV_def[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
  b <- mvrnorm(N, mu = A, Sigma = D)
  b <- b[rep(1:nrow(b), each = 8), ]
  y_hat <- Likelihoods_CV(fc = fc, b = b)
  prob[draw] <- mean(y_hat)
  hitrate[draw] <- sum(y_hat > 0.5)
  
} # This might take a few minutes

RespObsTable["CV_def", "HitRate"] <- mean(hitrate) / (N * 8)
RespObsTable["CV_def", "AvgP"]    <- mean(prob)

# Without Covariates (default settings)
hitrate <- rep(NA, nDraws)
prob <- rep(NA, nDraws)
fc <- NULL
A <- colMeans(NCV_def[["A"]][, -1])
D <- apply(X = NCV_def[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
  b <- mvrnorm(N, mu = A, Sigma = D)
  b <- b[rep(1:nrow(b), each = 8), ]
  y_hat <- Likelihoods_NCV(fc = fc, b = b)
  prob[draw] <- mean(y_hat)
  hitrate[draw] <- sum(y_hat > 0.5)
  
} # This might take a few minutes

RespObsTable["NCV_def", "HitRate"] <- mean(hitrate) / (N * 8)
RespObsTable["NCV_def", "AvgP"]    <- mean(prob)

# With Covariates (modified settings)
hitrate <- rep(NA, nDraws)
prob <- rep(NA, nDraws)
fc <- colMeans(CV_mod[["F"]][, -1])
A <- colMeans(CV_mod[["A"]][, -1])
D <- apply(X = CV_mod[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
  b <- mvrnorm(N, mu = A, Sigma = D)
  b <- b[rep(1:nrow(b), each = 8), ]
  y_hat <- Likelihoods_CV(fc = fc, b = b)
  prob[draw] <- mean(y_hat)
  hitrate[draw] <- sum(y_hat > 0.5)
  
} # This might take a few minutes

RespObsTable["CV_mod", "HitRate"] <- mean(hitrate) / (N * 8)
RespObsTable["CV_mod", "AvgP"]    <- mean(prob)

# Without Covariates (default settings)
hitrate <- rep(NA, nDraws)
prob <- rep(NA, nDraws)
fc <- NULL
A <- colMeans(NCV_mod[["A"]][, -1])
D <- apply(X = NCV_mod[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
  b <- mvrnorm(N, mu = A, Sigma = D)
  b <- b[rep(1:nrow(b), each = 8), ]
  y_hat <- Likelihoods_NCV(fc = fc, b = b)
  prob[draw] <- mean(y_hat)
  hitrate[draw] <- sum(y_hat > 0.5)
  
} # This might take a few minutes

RespObsTable["NCV_mod", "HitRate"] <- mean(hitrate) / (N * 8)
RespObsTable["NCV_mod", "AvgP"]    <- mean(prob)

##############################################
# Single Task Hold-Out Sample

choicedata <- cd[sequence == 3 & out_hold == 0]
SingleObsTable <- data.frame(HitRate = rep(NA, 4),
                             AvgP = rep(NA, 4),
                             row.names = c("CV_def", "NCV_def", "CV_mod", "NCV_mod"))

# With Covariates (default settings)
fc <- colMeans(CV_def[["F"]][, -1])
b <- CV_def[["C"]][, -c(1:2)]

y_hat <- Likelihoods_CV(fc = fc, b = b)
SingleObsTable["CV_def", "HitRate"] <- sum(y_hat > 0.5) / nrow(b)
SingleObsTable["CV_def", "AvgP"]    <- mean(y_hat)

# Without Covariates (default settings)
fc <- NULL
b <- NCV_def[["C"]][, -c(1:2)]

y_hat <- Likelihoods_NCV(fc = fc, b = b)
SingleObsTable["NCV_def", "HitRate"] <- sum(y_hat > 0.5) / nrow(b)
SingleObsTable["NCV_def", "AvgP"]    <- mean(y_hat)

# With Covariates (modified settings)
fc <- colMeans(CV_mod[["F"]][, -1])
b <- CV_mod[["C"]][, -c(1:2)]

y_hat <- Likelihoods_CV(fc = fc, b = b)
SingleObsTable["CV_mod", "HitRate"] <- sum(y_hat > 0.5) / nrow(b)
SingleObsTable["CV_mod", "AvgP"]    <- mean(y_hat)

# Without Covariates (modified settings)
fc <- NULL
b <- NCV_mod[["C"]][, -c(1:2)]

y_hat <- Likelihoods_NCV(fc = fc, b = b)
SingleObsTable["NCV_mod", "HitRate"] <- sum(y_hat > 0.5) / nrow(b)
SingleObsTable["NCV_mod", "AvgP"]    <- mean(y_hat)

##############################################
# Result Tables

RLH_table
LLSimTable
RespObsTable
SingleObsTable
