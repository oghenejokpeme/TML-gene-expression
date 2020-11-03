#!/usr/bin/Rscript --vanilla
source("lib_general_functions.R")
source("lib_general_learners.R")

PerformRoundA <- function(all.genes, parallel_ = TRUE){ 
  for (gene in all.genes) { 
    print(paste0("Building for: ", gene))
    # Loag gene base data.
    gene.df <- LoadBaseGeneData(gene)
    train <- as.h2o(gene.df$train)
    test <- as.h2o(gene.df$test)
    
    # Train and log model. 
    model.id <- paste0(gene, "_nn_RA")
    nn.model <- FitNeuralNet(train, model.id)
    
    save.path <- "../output/1_transformative_learning/models_ra/" 
    h2o.saveModel(object = nn.model, path = save.path, force = TRUE) 
    
    # Get performance on test set.
    # dl.perf <- h2o.performance(model = nn.fit, newdata = test)      
    test.predictions <- h2o.predict(nn.model, test)
    test.predictions <- as.vector(test.predictions[, 1])
    nn.perf <- GetPerformanceMetrics(gene.df$test[, c("y")], test.predictions)

    # Write performance metrics to file.
    perf.metrics <- c(nn.perf$rsquared, nn.perf$mse, nn.perf$rmse)
    gene.perf <- c(gene, perf.metrics) 
    write(gene.perf, 
          ncolumns = length(gene.perf),
          append = TRUE, 
          file = paste0("../output/1_transformative_learning/",
                        "results/base_nn_results.txt")
          )
  }
}

PerformExperiments <- function(parallel_ = FALSE){
  # Setup H2O.
  h2o.init(nthreads = -1) 
  h2o.no_progress()
  # Load gene list.
  all.genes <- LoadGeneList()
  PerformRoundA(all.genes) 
}

PerformExperiments()
