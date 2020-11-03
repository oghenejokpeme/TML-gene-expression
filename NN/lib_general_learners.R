#!/usr/bin/Rscript --vanilla

library(h2o)
library(doMC)

FitNeuralNet <- function(train, model.id){
  cap <- dim(train)[2] - 1
  x <- colnames(train)[1:cap]
  y <- "y"
  
  flayer <- ceiling((1/3) * length(x))
  slayer <- ceiling((1/3) * flayer)
    
  nn.fit <- h2o.deeplearning(x = x,
                             y = y,
                             model_id = model.id, 
                             training_frame = train,
                             standardize = FALSE,
                             activation = "Tanh",
                             hidden = c(flayer, slayer),
                             seed = 35542)

  return(nn.fit)
}