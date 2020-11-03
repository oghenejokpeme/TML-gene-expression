#!/usr/bin/Rscript --vanilla

library(parallel)
library(doMC)
library(doRNG)
library(foreach)
library(e1071)

FitSvmReg <- function(x, y, parallel_=FALSE){
  ctrl <- trainControl(method = "cv", number = 5, allowParallel = parallel_)
  svmgrid <- expand.grid(sigma = 2^seq(-15, 0, by = 1),
                         C = c(0.1, 1, 10))
  set.seed(53475)
  svm.selection <- train(x = x, 
                         y = y,
                         scaled = FALSE,
                         method = "svmRadial",
                         metric = "RMSE", 
                         epsilon = 0.01,
                         trControl = ctrl,
                         tuneGrid = svmgrid) 
  svm.model <- svm.selection$finalModel
    
  return(svm.model)
}

FitSingleSVM <- function(x, y){
  svm.model <- svm(x = x, y = y, scale = FALSE, 
                   type = "eps-regression", epsilon = 0.01,
                   cost = 0.25, gamma = 0.5)
  
  return(svm.model)
}