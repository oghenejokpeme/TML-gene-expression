#!/usr/bin/Rscript --vanilla
source("lib_general_functions.R")
source("lib_general_learners.R")

PerformRoundA <- function(all.genes, parallel_ = TRUE){
  for (gene in all.genes) {
    # Loag gene base data.
    gene.df <- LoadBaseGeneData(gene)
    x.train <- gene.df$x.train
    y.train <- gene.df$y.train
    x.test <- gene.df$x.test
    y.test <- gene.df$y.test

    # Train and log model.
    rf.model <- FitRandomForest(x.train, y.train)
    saveRDS(rf.model,
            file = paste0("../output/1_transformative_learning/",
                          "models_ra/", gene, "_rf_RA.rds"))

    # Get performance on test set.
    test.predictions <- predict(rf.model, x.test)
    rf.perf <- GetPerformanceMetrics(y.test, test.predictions)

    # Write performance metrics to file.
    perf.metrics <- c(rf.perf$rsquared, rf.perf$mse, rf.perf$rmse)
    gene.perf <- c(gene, perf.metrics) 
    write(gene.perf, 
          ncolumns = length(gene.perf),
          append = TRUE, 
          file = paste0("../output/1_transformative_learning/",
                        "results/base_rf_results.txt"))
  }
}

PerformExperiments <- function(parallel_ = TRUE){
  if (parallel_){ 
    registerDoMC(cores = 50) 
  }
  # Load gene list.
  all.genes <- LoadGeneList()
  PerformRoundA(all.genes) 
}

PerformExperiments()
