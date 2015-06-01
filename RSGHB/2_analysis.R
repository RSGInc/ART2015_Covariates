##############################################
# Load packages, data, and model objects

setwd("C://Users/jeff.keller/Desktop/ART2015_Covariates/RSGHB/")

library(RSGHB)
library(data.table)
library(MASS)

source("0_likelihood_functions.R")

# Load RSGHB results
load("CV_def.RData")
load("CV_mod.RData")
load("NCV_def.RData")
load("NCV_mod.RData")
load("choicedata.RData")

# Load CBC HB results
NCV_saw <- list()
CV_saw <- list()
NCV_saw$A <- fread("../CBCHB/1_NCV_alpha.csv")
NCV_saw$C <- fread("../CBCHB/1_NCV_utilities.csv")
NCV_saw$D <- fread("../CBCHB/1_NCV_covariances.csv")
CV_saw$A <- fread("../CBCHB/2_CV_alpha.csv", skip = 2) # has an extra header row for the covariate descriptions
CV_saw$C <- fread("../CBCHB/2_CV_utilities.csv")
CV_saw$D <- fread("../CBCHB/2_CV_covariances.csv")

##############################################
# Organize CBC HB results

# Drop the burn-in iterations
NCV_saw$A <- NCV_saw$A[-c(1:100000)]
CV_saw$A  <-  CV_saw$A[-c(1:100000)]
NCV_saw$D <- NCV_saw$D[-c(1:100000)]
CV_saw$D  <-  CV_saw$D[-c(1:100000)]
gc()

# Convert to matrices
NCV_saw$A <- as.matrix(NCV_saw$A)
NCV_saw$C <- as.matrix(NCV_saw$C)
NCV_saw$D <- as.matrix(NCV_saw$D)[, -1]
CV_saw$A  <- as.matrix(CV_saw$A)
CV_saw$C  <- as.matrix(CV_saw$C)
CV_saw$D  <- as.matrix(CV_saw$D)[, -1]

# Drop the baseline columns
NCV_saw$A <- NCV_saw$A[, apply(X = NCV_saw$A, MARGIN = 2, FUN = function(x) !all(x == 0))]
NCV_saw$C <- NCV_saw$C[, apply(X = NCV_saw$C, MARGIN = 2, FUN = function(x) !all(x == 0))]
CV_saw$A  <-  CV_saw$A[, apply(X =  CV_saw$A, MARGIN = 2, FUN = function(x) !all(x == 0))]
CV_saw$C  <-  CV_saw$C[, apply(X =  CV_saw$C, MARGIN = 2, FUN = function(x) !all(x == 0))]

# Convert covariance matrices into 3D arrays
tmp <- array(dim = c(27, 27, 10000))
for (i in 1:nrow(NCV_saw$D)) {
  mat <- matrix(NA, nrow = 27, ncol = 27)
  mat[lower.tri(mat, diag = TRUE)] <- NCV_saw$D[i, ]
  mat[upper.tri(mat)] <- t(mat)[upper.tri(t(mat))]
  tmp[,, i] <- mat
}
NCV_saw$D <- tmp

tmp <- array(dim = c(27, 27, 10000))
for (i in 1:nrow(CV_saw$D)) {
  mat <- matrix(NA, nrow = 27, ncol = 27)
  mat[lower.tri(mat, diag = TRUE)] <- CV_saw$D[i, ]
  mat[upper.tri(mat)] <- t(mat)[upper.tri(t(mat))]
  tmp[,, i] <- mat
}
CV_saw$D <- tmp
rm(tmp)
gc()

##############################################
# Model Summaries

# Lower-level model RLH
RLH_table <- data.frame(RLH = c(mean(CV_def [["C"]][, "RLH"]),
                                mean(NCV_def[["C"]][, "RLH"]),
                                mean(CV_mod [["C"]][, "RLH"]),
                                mean(NCV_mod[["C"]][, "RLH"]),
                                mean(CV_saw [["C"]][, "RLH"]),
                                mean(NCV_saw[["C"]][, "RLH"])),
                        row.names = c("CV_def", "NCV_def", "CV_mod", "NCV_mod", "CV_saw", "NCV_saw"))

##############################################
# Upper-level model simulated LL

N <- 404
nDraws <- 1000
choicedata <- cd[out_hold == 0 & sequence != 3]
LLSimTable <- data.frame(AvgLL = rep(NA, 6), row.names = c("CV_def", "NCV_def", "CV_mod", "NCV_mod", "CV_saw", "NCV_saw"))

# RSGHB With Covariates (default settings)
ll <- matrix(NA, N, nDraws)
fc <- colMeans(CV_def[["F"]][, -1])
A <- colMeans(CV_def[["A"]][, -1])
D <- apply(X = CV_def[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
     b <- mvrnorm(N, mu = A, Sigma = D)
     b <- b[rep(1:nrow(b), each = 7), ]
     ll[, draw] <- aggregate(Likelihoods_CV(fc = fc, b = b), list(rep(1:N, each = 7)), prod)[, 2]
     if ((draw %% 25) == 0) print(draw)
} # This might take a few minutes

LLSimTable["CV_def", "AvgLL"] <- sum(log(rowMeans(ll)))

# RSGHB Without Covariates (default settings)
ll <- matrix(NA, N, nDraws)
fc <- NULL
A <- colMeans(NCV_def[["A"]][, -1])
D <- apply(X = NCV_def[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
     b <- mvrnorm(N, mu = A, Sigma = D)
     b <- b[rep(1:nrow(b), each = 7), ]
     ll[, draw] <- aggregate(Likelihoods_NCV(fc = fc, b = b), list(rep(1:N, each = 7)), prod)[, 2]
     if ((draw %% 25) == 0) print(draw)
} # This might take a few minutes

LLSimTable["NCV_def", "AvgLL"] <- sum(log(rowMeans(ll)))

# RSGHB With Covariates (modified settings)
ll <- matrix(NA, N, nDraws)
fc <- colMeans(CV_mod[["F"]][, -1])
A <- colMeans(CV_mod[["A"]][, -1])
D <- apply(X = CV_mod[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
     b <- mvrnorm(N, mu = A, Sigma = D)
     b <- b[rep(1:nrow(b), each = 7), ]
     ll[, draw] <- aggregate(Likelihoods_CV(fc = fc, b = b), list(rep(1:N, each = 7)), prod)[, 2]
     if ((draw %% 25) == 0) print(draw)
} # This might take a few minutes

LLSimTable["CV_mod", "AvgLL"] <- sum(log(rowMeans(ll)))

# RSGHB Without Covariates (modified settings)
ll <- matrix(NA, N, nDraws)
fc <- NULL
A <- colMeans(NCV_mod[["A"]][, -1])
D <- apply(X = NCV_mod[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
     b <- mvrnorm(N, mu = A, Sigma = D)
     b <- b[rep(1:nrow(b), each = 7), ]
     ll[, draw] <- aggregate(Likelihoods_NCV(fc = fc, b = b), list(rep(1:N, each = 7)), prod)[, 2]
     if ((draw %% 25) == 0) print(draw)
} # This might take a few minutes

LLSimTable["NCV_mod", "AvgLL"] <- sum(log(rowMeans(ll)))

# CBC HB With Covariates
ll <- matrix(NA, N, nDraws)
fc <- colMeans(CV_saw[["A"]][, rep(29:55, each = 2) + rep(c(0, 27), times = 27)]) # CBC HB stores covariates with the means and in a different order
A <- colMeans(CV_saw[["A"]][, 2:28])
D <- apply(X = CV_saw[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
  b <- mvrnorm(N, mu = A, Sigma = D)
  b <- b[rep(1:nrow(b), each = 7), ]
  ll[, draw] <- aggregate(Likelihoods_CV(fc = fc, b = b), list(rep(1:N, each = 7)), prod)[, 2]
  if ((draw %% 25) == 0) print(draw)
} # This might take a few minutes

LLSimTable["CV_saw", "AvgLL"] <- sum(log(rowMeans(ll)))

# CBC HB Without Covariates
ll <- matrix(NA, N, nDraws)
fc <- NULL
A <- colMeans(NCV_saw[["A"]][, -1])
D <- apply(X = NCV_saw[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
  b <- mvrnorm(N, mu = A, Sigma = D)
  b <- b[rep(1:nrow(b), each = 7), ]
  ll[, draw] <- aggregate(Likelihoods_NCV(fc = fc, b = b), list(rep(1:N, each = 7)), prod)[, 2]
  if ((draw %% 25) == 0) print(draw)
} # This might take a few minutes

LLSimTable["NCV_saw", "AvgLL"] <- sum(log(rowMeans(ll)))

##############################################
# Complete Respondent Hold-Out Sample

N <- 45
nDraws <- 1000
choicedata <- cd[out_hold == 1]
RespObsTable <- data.frame(HitRate = rep(NA, 6),
                           AvgP = rep(NA, 6),
                           LL   = rep(NA, 6),
                           row.names = c("CV_def", "NCV_def", "CV_mod", "NCV_mod", "CV_saw", "NCV_saw"))

# RSGHB With Covariates (default settings)
hitrate <- matrix(NA, N*8, nDraws)
prob <- matrix(NA, N*8, nDraws)
ll <- matrix(NA, N, nDraws)
fc <- colMeans(CV_def[["F"]][, -1])
A <- colMeans(CV_def[["A"]][, -1])
D <- apply(X = CV_def[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
     b <- mvrnorm(N, mu = A, Sigma = D)
     b <- b[rep(1:nrow(b), each = 8), ]
     y_hat <- Likelihoods_CV(fc = fc, b = b)
     ll[, draw] <- aggregate(y_hat, list(rep(1:N, each = 8)), prod)[, 2]     
     prob[, draw] <- y_hat
     hitrate[, draw] <- y_hat > 0.5
     if ((draw %% 25) == 0) print(draw)     
} # This might take a few minutes

RespObsTable["CV_def", "HitRate"] <- mean(hitrate) 
RespObsTable["CV_def", "AvgP"]    <- mean(prob)
RespObsTable["CV_def", "LL"]    <- sum(log(rowMeans(ll)))

# RSGHB Without Covariates (default settings)
hitrate <- matrix(NA, N*8, nDraws)
prob <- matrix(NA, N*8, nDraws)
ll <- matrix(NA, N, nDraws)
fc <- NULL
A <- colMeans(NCV_def[["A"]][, -1])
D <- apply(X = NCV_def[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
     b <- mvrnorm(N, mu = A, Sigma = D)
     b <- b[rep(1:nrow(b), each = 8), ]
     y_hat <- Likelihoods_NCV(fc = fc, b = b)
     ll[, draw] <- aggregate(y_hat, list(rep(1:N, each = 8)), prod)[, 2]     
     prob[, draw] <- y_hat
     hitrate[, draw] <- y_hat > 0.5
     if ((draw %% 25) == 0) print(draw)  
} # This might take a few minutes

RespObsTable["NCV_def", "HitRate"] <- mean(hitrate) 
RespObsTable["NCV_def", "AvgP"]    <- mean(prob)
RespObsTable["NCV_def", "LL"]      <- sum(log(rowMeans(ll)))

# RSGHB With Covariates (modified settings)
hitrate <- matrix(NA, N*8, nDraws)
prob <- matrix(NA, N*8, nDraws)
ll <- matrix(NA, N, nDraws)
fc <- colMeans(CV_mod[["F"]][, -1])
A <- colMeans(CV_mod[["A"]][, -1])
D <- apply(X = CV_mod[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
     b <- mvrnorm(N, mu = A, Sigma = D)
     b <- b[rep(1:nrow(b), each = 8), ]
     y_hat <- Likelihoods_CV(fc = fc, b = b)
     ll[, draw] <- aggregate(y_hat, list(rep(1:N, each = 8)), prod)[, 2]     
     prob[, draw] <- y_hat
     hitrate[, draw] <- y_hat > 0.5
     if ((draw %% 25) == 0) print(draw)     
} # This might take a few minutes

RespObsTable["CV_mod", "HitRate"] <- mean(hitrate) 
RespObsTable["CV_mod", "AvgP"]    <- mean(prob)
RespObsTable["CV_mod", "LL"]      <- sum(log(rowMeans(ll)))

# RSGHB Without Covariates (default settings)
hitrate <- matrix(NA, N*8, nDraws)
prob <- matrix(NA, N*8, nDraws)
ll <- matrix(NA, N, nDraws)
A <- colMeans(NCV_mod[["A"]][, -1])
D <- apply(X = NCV_mod[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
     b <- mvrnorm(N, mu = A, Sigma = D)
     b <- b[rep(1:nrow(b), each = 8), ]
     y_hat <- Likelihoods_NCV(fc = fc, b = b)
     ll[, draw] <- aggregate(y_hat, list(rep(1:N, each = 8)), prod)[, 2]     
     prob[, draw] <- y_hat
     hitrate[, draw] <- y_hat > 0.5
     if ((draw %% 25) == 0) print(draw)     
} # This might take a few minutes

RespObsTable["NCV_mod", "HitRate"] <- mean(hitrate) 
RespObsTable["NCV_mod", "AvgP"]    <- mean(prob)
RespObsTable["NCV_mod", "LL"]      <- sum(log(rowMeans(ll)))

# CBC HB With Covariates
hitrate <- matrix(NA, N*8, nDraws)
prob <- matrix(NA, N*8, nDraws)
ll <- matrix(NA, N, nDraws)
fc <- colMeans(CV_saw[["A"]][, rep(29:55, each = 2) + rep(c(0, 27), times = 27)]) # CBC HB stores covariates with the means and in a different order
A <- colMeans(CV_saw[["A"]][, 2:28])
D <- apply(X = CV_saw[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
  b <- mvrnorm(N, mu = A, Sigma = D)
  b <- b[rep(1:nrow(b), each = 8), ]
  y_hat <- Likelihoods_CV(fc = fc, b = b)
  ll[, draw] <- aggregate(y_hat, list(rep(1:N, each = 8)), prod)[, 2]     
  prob[, draw] <- y_hat
  hitrate[, draw] <- y_hat > 0.5
  if ((draw %% 25) == 0) print(draw)     
} # This might take a few minutes

RespObsTable["CV_saw", "HitRate"] <- mean(hitrate) 
RespObsTable["CV_saw", "AvgP"]    <- mean(prob)
RespObsTable["CV_saw", "LL"]      <- sum(log(rowMeans(ll)))

# CBC HB Without Covariates
hitrate <- matrix(NA, N*8, nDraws)
prob <- matrix(NA, N*8, nDraws)
ll <- matrix(NA, N, nDraws)
fc <- NULL
A <- colMeans(CV_saw[["A"]][, 2:28])
D <- apply(X = CV_saw[["D"]], MARGIN = 1:2, FUN = mean)

set.seed(1987)
for (draw in 1:nDraws) {
  b <- mvrnorm(N, mu = A, Sigma = D)
  b <- b[rep(1:nrow(b), each = 8), ]
  y_hat <- Likelihoods_NCV(fc = fc, b = b)
  ll[, draw] <- aggregate(y_hat, list(rep(1:N, each = 8)), prod)[, 2]     
  prob[, draw] <- y_hat
  hitrate[, draw] <- y_hat > 0.5
  if ((draw %% 25) == 0) print(draw)     
} # This might take a few minutes

RespObsTable["NCV_saw", "HitRate"] <- mean(hitrate) 
RespObsTable["NCV_saw", "AvgP"]    <- mean(prob)
RespObsTable["NCV_saw", "LL"]      <- sum(log(rowMeans(ll)))

##############################################
# Single Task Hold-Out Sample

choicedata <- cd[sequence == 3 & out_hold == 0]
SingleObsTable <- data.frame(HitRate = rep(NA, 6),
                             AvgP = rep(NA, 6),
                             row.names = c("CV_def", "NCV_def", "CV_mod", "NCV_mod", "CV_saw", "NCV_saw"))

# RSGHB With Covariates (default settings)
fc <- colMeans(CV_def[["F"]][, -1])
b <- CV_def[["C"]][, -c(1:2)]

y_hat <- Likelihoods_CV(fc = fc, b = b)
SingleObsTable["CV_def", "HitRate"] <- sum(y_hat > 0.5) / nrow(b)
SingleObsTable["CV_def", "AvgP"]    <- mean(y_hat)

# RSGHB Without Covariates (default settings)
fc <- NULL
b <- NCV_def[["C"]][, -c(1:2)]

y_hat <- Likelihoods_NCV(fc = fc, b = b)
SingleObsTable["NCV_def", "HitRate"] <- sum(y_hat > 0.5) / nrow(b)
SingleObsTable["NCV_def", "AvgP"]    <- mean(y_hat)

# RSGHB With Covariates (modified settings)
fc <- colMeans(CV_mod[["F"]][, -1])
b <- CV_mod[["C"]][, -c(1:2)]

y_hat <- Likelihoods_CV(fc = fc, b = b)
SingleObsTable["CV_mod", "HitRate"] <- sum(y_hat > 0.5) / nrow(b)
SingleObsTable["CV_mod", "AvgP"]    <- mean(y_hat)

# RSGHB Without Covariates (modified settings)
fc <- NULL
b <- NCV_mod[["C"]][, -c(1:2)]

y_hat <- Likelihoods_NCV(fc = fc, b = b)
SingleObsTable["NCV_mod", "HitRate"] <- sum(y_hat > 0.5) / nrow(b)
SingleObsTable["NCV_mod", "AvgP"]    <- mean(y_hat)

# CBC HB With Covariates
fc <- NULL # CBC HB bakes the covariates into the resulting individual utilities
b <- CV_saw[["C"]][, -c(1:2)]

y_hat <- Likelihoods_NCV(fc = fc, b = b)
SingleObsTable["CV_saw", "HitRate"] <- sum(y_hat > 0.5) / nrow(b)
SingleObsTable["CV_saw", "AvgP"]    <- mean(y_hat)

# CBC HB Without Covariates
fc <- NULL
b <- NCV_saw[["C"]][, -c(1:2)]

y_hat <- Likelihoods_NCV(fc = fc, b = b)
SingleObsTable["NCV_saw", "HitRate"] <- sum(y_hat > 0.5) / nrow(b)
SingleObsTable["NCV_saw", "AvgP"]    <- mean(y_hat)

##############################################
# Result Tables

RLH_table
LLSimTable
RespObsTable
SingleObsTable


# ADD KS TESTS between def/mod models and secondarily between RSG vs. CBCHB