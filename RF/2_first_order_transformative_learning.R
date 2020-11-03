#!/usr/bin/Rscript --vanilla
source("lib_general_functions.R")
source("lib_general_learners.R")

WriteChunksToFile <- function(gene, chunked.metagenes){
  base.path <- paste0("../output/1_transformative_learning/",
                      "chunks/", gene)
  dir.create(base.path, showWarnings = FALSE)

  for (i in 1:length(chunked.metagenes)) {
    chunk <- chunked.metagenes[[i]]
    file.path <- paste0(base.path, "/", i, ".txt")
    write(chunk, file= file.path, ncolumns = 1)
  }
}

# nin.chunks is number of meta-gene predictions in each chunk.
CreateMetaGeneList <- function(all.genes, nin.chunks = 50){
  genes.metagenes <- NULL
  no.chunks <- NULL
  for (gene in all.genes) {
    gene.metagenes <- setdiff(all.genes, gene)
    x <- seq_along(gene.metagenes)
    chunked.metagenes <- split(gene.metagenes, ceiling(x/nin.chunks))
    genes.metagenes[[gene]] <- chunked.metagenes

    if (is.null(no.chunks)) {
      no.chunks <- length(chunked.metagenes)
    }

    WriteChunksToFile(gene, chunked.metagenes)
  }
  output <- list(no.chunks = no.chunks, genes.metagenes = genes.metagenes)

  return(output)
}

CreateGeneLoadOrder <- function(all.genes, nin.chunks = 50){
  x <- seq_along(all.genes)
  temp.loadorder <- split(all.genes, ceiling(x/nin.chunks))
  total.chunks <- length(temp.loadorder)

  genes.loadorder <- NULL
  for (i in 1:total.chunks) {
    if (i == 1) {
      genes.loadorder[[i]] <- c(temp.loadorder[[i]], temp.loadorder[[i + 1]][1]) 
    } else if ( i == total.chunks) {
      genes.loadorder[[i]] <- c(temp.loadorder[[i - 1]][nin.chunks], 
                                temp.loadorder[[i]])
    } else {
      genes.loadorder[[i]] <- c(temp.loadorder[[i - 1]][nin.chunks], 
                                temp.loadorder[[i]], 
                                temp.loadorder[[i + 1]][1]
                                )
    }
  }

  return(genes.loadorder)
}

CopyModelsToMemory <- function(all.genes){
  gene.models <- list()
  print("Copying models to memory ...") 
  for (gene in all.genes) { 
    print(paste0(" Load ", which(all.genes == gene), " ", gene))
    model.path <- paste0("../output/1_transformative_learning/",
                         "models_ra/", gene, "_rf_RA.rds")
    model <- readRDS(model.path)
    gene.models[[gene]] <- model
  }

  return(gene.models)
}

PredictNewInputs <- function(cgene, gene.df, current.metagenes, gene.models, 
                             chunk.no, parallel_ = TRUE){
  print(paste0("Predicting news descriptors for ", cgene, " chunk ", chunk.no))
  transx.train <- NULL
  transx.test <- NULL

  if (parallel_) {
    transx.train <- foreach (gene =  current.metagenes, .combine = 'cbind', 
                             .inorder = TRUE) %dopar% {
                      print(paste0("  +", 
                                   which(current.metagenes == gene), 
                                   " ", 
                                   gene
                                  )
                           )
                      temp.train <- predict(gene.models[[gene]], gene.df$x.train)
                    }
    transx.test <- foreach (gene =  current.metagenes, .combine = 'cbind', 
                            .inorder = TRUE) %dopar% {
                      print(paste0("  -", 
                                   which(current.metagenes == gene), 
                                   " ", 
                                   gene
                                  )
                          )
                      temp.test <- predict(gene.models[[gene]], gene.df$x.test)
                   } 
  } else {
    for (gene in current.metagenes) {      
      print(paste0("  ", which(current.metagenes == gene), " ", gene))
      
      temp.train <- predict(gene.models[[gene]], gene.df$x.train)
      temp.test <- predict(gene.models[[gene]], gene.df$x.test)
      transx.train <- cbind(transx.train, temp.train)
      transx.test <- cbind(transx.test, temp.test)
    }
  }
  
  colnames(transx.train) <- current.metagenes
  colnames(transx.test) <- current.metagenes

  # Save transformed input to file.
  write.csv(transx.train,
            file = paste0("../output/1_transformative_learning/",
                          "transformed_inputs/", cgene, 
                          "_", chunk.no, "_rf_first_order_train.csv")
            )
  write.csv(transx.test,
            file = paste0("../output/1_transformative_learning/",
                          "transformed_inputs/", cgene, 
                          "_", chunk.no, "_rf_first_order_test.csv")
            )
}

LogExperimentFeatures <- function(gene, chunk.no, df.x){
  base.path <- paste0("../output/1_transformative_learning/",
                      "build_logs/", gene, "_rf.txt")

  df.dim <- dim(df.x$train)
  train.metagenes <- colnames(df.x$train)
  test.metagenes <- colnames(df.x$test)


  write(c(chunk.no, df.dim), 
        file = base.path, 
        append = TRUE, 
        ncolumns = length(df.dim) + 1
       )
  write(train.metagenes, 
        file = base.path, 
        append = TRUE, 
        ncolumns = length(train.metagenes)
       )
  write(test.metagenes, 
        file = base.path, 
        append = TRUE, 
        ncolumns = length(test.metagenes)
       )
  write("", file = base.path, append = TRUE, ncolumns = 1)
}

PerformRoundB <- function(all.genes, metaplan.out, genes.loadorder, 
                          parallel_ = TRUE){
  no.chunks <- metaplan.out$no.chunks
  genes.metagenes <- metaplan.out$genes.metagenes
 
  for (chunk.no in 1:10) { 
    print("")
    print(paste0("Performing chunk ", chunk.no, " experiments."))
    # Load models required for this chunk to memory.
    chunk.models <- genes.loadorder[[chunk.no]]
    gene.models <- CopyModelsToMemory(chunk.models)
    
    for (gene in all.genes) { 
      print(paste0("*", gene))
      gene.df <- LoadBaseGeneData(gene)
      
      # Predict new features.
      current.metagenes <- genes.metagenes[[gene]][[chunk.no]]
      PredictNewInputs(gene, gene.df, current.metagenes, gene.models, chunk.no)
      
      # Load latest data for predicted features.
      df.x <- LoadFirstOrderInputs(gene, chunk.no)
      LogExperimentFeatures(gene, chunk.no, df.x)

      y.train <- gene.df$y.train
      y.test <- gene.df$y.test

      # Build and save new models and results.
      rf.model <- FitRandomForest(df.x$train, y.train)

      # Get performance on test set.
      test.predictions <- predict(rf.model, df.x$test)
      rf.perf <- GetPerformanceMetrics(y.test, test.predictions)

      # Write performance metrics to file.
      perf.metrics <- c(rf.perf$rsquared, rf.perf$mse, rf.perf$rmse)
      gene.perf <- c(gene, perf.metrics) 
      
      saveRDS(rf.model,
              file = paste0("../output/1_transformative_learning/",
                            "models_rb/", gene, "_", chunk.no, 
                            "_rf_RB.rds")
              )

      write(gene.perf, 
            ncolumns = length(gene.perf),
            append = TRUE, 
            file = paste0("../output/1_transformative_learning/",
                          "results/firstorder_rf_", chunk.no, 
                          "_results.txt")
            )
    }
  }
}

PerformExperiments <- function(parallel_ = TRUE, nin.chunks = 50){
  if (parallel_){ 
    registerDoMC(cores = 50) 
  }
  # Load gene list.
  all.genes <- LoadGeneList()

  genes.loadorder <- CreateGeneLoadOrder(all.genes, nin.chunks)
  metaplan.out <- CreateMetaGeneList(all.genes, nin.chunks)
  PerformRoundB(all.genes, metaplan.out, genes.loadorder) 
}

PerformExperiments()
