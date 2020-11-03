#!/usr/bin/Rscript --vanilla

library(data.table)
library(parallel)
library(doMC)
library(foreach)

# Get general peformance metrics. 
GetPerformanceMetrics <- function(y.actual, y.predicted){
  rsquared <- 1 - (sum((y.actual - y.predicted)^2) / 
                   sum((y.actual - mean(y.actual))^2))
  mse <- sum((y.actual - y.predicted)^2) / length(y.actual)
  rmse <- sqrt(sum((y.actual - y.predicted)^2) / length(y.actual))
    
  performance <- list(rsquared = rsquared,
                      mse = mse,
                      rmse = rmse)

  return(performance)
}

LoadGeneList <- function(){
  all.genes <- read.table(file = ".../input/data_buckets/all.genes.txt",
                          header = FALSE)[,1]
  all.genes <- as.vector(all.genes)

  return(all.genes)
}

LoadBaseGeneData <- function(gene){
  base.path <- paste0("../input/experiment_data/", gene, "_")
  df.train <- read.csv(file = paste0(base.path, "base_train.csv"),
                       header = TRUE,
                       row.names = 1)
  df.test <- read.csv(file = paste0(base.path, "base_test.csv"),
                      header = TRUE,
                      row.names = 1)

  df <- list(train = df.train,
             test = df.test)

  return(df)
}

LoadFirstOrderInputs <- function(gene, chunk.no){
  df.train <- NULL
  df.test <- NULL

  for (i in 1:chunk.no) {
    tf.train <- data.matrix(data.frame(fread(paste0("../output/1_transformative_learning/",
                                                  "transformed_inputs/", gene, 
                                                  "_", i, "_nn_first_order_train.csv"), 
                                           header = TRUE),
                                    row.names = 1))
    tf.test <- data.matrix(data.frame(fread(paste0("../output/1_transformative_learning/",
                                                 "transformed_inputs/", gene, 
                                                  "_", i, "_nn_first_order_test.csv"), 
                                          header = TRUE),
                                    row.names = 1)) 
    df.train <- cbind(df.train, tf.train)
    df.test <- cbind(df.test, tf.test)
  }
  
  df <- list(train = df.train, test = df.test)

  return(df)
}
