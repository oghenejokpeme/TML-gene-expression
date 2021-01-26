library(caret)
library(doMC)
library(doRNG)
library(foreach)
library(parallel)
library(xgboost)

FitXgb <- function(x, y){
  ctrl <- trainControl(
    method = "LGOCV", 
    p = 0.7,
    number = 1,
    allowParallel = FALSE
  )

  grid <- expand.grid(
    nrounds = c(1000, 1500),
    max_depth = 6,
    eta = c(0.001, 0.01, 0.1, 0.2, 0.3),
    gamma = 0.05,
    subsample = 0.5,
    colsample_bytree = 1,
    min_child_weight = 1
  )

  # Having objective = "reg:squarederror" below causes the following
  # warning without any detriment to actual model building. 
  # "The following parameters were provided multiple times:
	#    objective
  #  Only the last value for each of them will be used."
  # If objective = "reg:squarederror" is not passed to train, the 
  # following warning is received:
  # "reg:linear is now deprecated in favor of reg:squarederror."
  # Therefore, the only option is to suppress local warnings within 
  # this function. - OIO, 31-Oct-2020 

  set.seed(47567)
  suppressWarnings(
    xgb.fit <- train(
      x = x, 
      y = y,
      method = "xgbTree",
      eval_metric = "rmse",
      verbose = 0, 
      trControl = ctrl,
      tuneGrid = grid,
      objective = "reg:squarederror"
    )
  )
  xgb.fit
}