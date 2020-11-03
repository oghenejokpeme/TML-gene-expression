#!/usr/bin/Rscript --vanilla
source("lib_general_functions.R")
source("lib_general_learners.R")

WriteChunksToFile <- function(gene, chunked.metagenes){
  base.path <- paste0("../output/1_transformative_learning/",
                      "chunks/", gene)
  dir.create(base.path, showWarnings = FALSE)

  for (i in 1:length(chunked.metagenes)) {
    chunk <- chunked.metagenes[[i]]
    file.path <- paste0(base.path, "/", i, "_nn.txt")
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
                         "models_ra/", gene, "_nn_RA")
    model <- h2o.loadModel(model.path)
    gene.models[[gene]] <- model
  }

  return(gene.models)
}

PredictNewInputs <- function(cgene, gene.df, current.metagenes, gene.models, 
                             chunk.no, parallel_ = FALSE){
  print(paste0("Predicting new descriptors for ", cgene, " chunk ", chunk.no))
  transx.train <- NULL
  transx.test <- NULL

  train <- as.h2o(gene.df$train)
  test <- as.h2o(gene.df$test)

  if (parallel_) {
    transx.train <- foreach (gene =  current.metagenes, .combine = 'cbind', 
                             .inorder = TRUE) %dopar% {
                      print(paste0("  +", 
                                   which(current.metagenes == gene), 
                                   " ", 
                                   gene
                                  )
                           )
                      temp.train <- as.vector(h2o.predict(gene.models[[gene]], 
                                                          train)[,1]
                                             )
                    }
    transx.test <- foreach (gene =  current.metagenes, .combine = 'cbind', 
                            .inorder = TRUE) %dopar% {
                      print(paste0("  -", 
                                   which(current.metagenes == gene), 
                                   " ", 
                                   gene
                                  )
                          )
                      temp.test <- as.vector(h2o.predict(gene.models[[gene]], 
                                                         test)[,1]
                                            )
                   } 
  } else {
    for (gene in current.metagenes) {      
      print(paste0("  ", which(current.metagenes == gene), " ", gene))
      
      temp.train <- as.vector(h2o.predict(gene.models[[gene]], train)[,1])
      temp.test <- as.vector(h2o.predict(gene.models[[gene]], test)[,1])
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
                          "_", chunk.no, "_nn_first_order_train.csv")
            )
  write.csv(transx.test,
            file = paste0("../output/1_transformative_learning/",
                          "transformed_inputs/", cgene, 
                          "_", chunk.no, "_nn_first_order_test.csv")
            )
}

LogExperimentFeatures <- function(gene, chunk.no, df.x){
  base.path <- paste0("../output/1_transformative_learning/",
                      "build_logs/", gene, "_nn.txt")

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
 
      y.train <- gene.df$train[, c("y")]
      y.test <- gene.df$test[, c("y")]  
       
      train <- as.h2o(cbind(df.x$train, y = y.train))
      test <- as.h2o(cbind(df.x$test, y = y.test))
            
      # Build and save new models and results.
      model.id <- paste0(gene, "_", chunk.no, "_nn_RB") 
      nn.model <- FitNeuralNet(train, model.id)
        
      save.path <- "../output/1_transformative_learning/models_rb/" 
      h2o.saveModel(object = nn.model, path = save.path, force = TRUE)
         
      # Get performance on test set.
      test.predictions <- h2o.predict(nn.model, test)
      test.predictions <- as.vector(test.predictions[,1])
      nn.perf <- GetPerformanceMetrics(y.test, test.predictions)

      # Write performance metrics to file.
      perf.metrics <- c(nn.perf$rsquared, nn.perf$mse, nn.perf$rmse)
      gene.perf <- c(gene, perf.metrics) 
             
      write(gene.perf, 
            ncolumns = length(gene.perf),
            append = TRUE, 
            file = paste0("../output/1_transformative_learning/",
                          "results/firstorder_nn_", chunk.no, 
                          "_results.txt")
            )
    }
  }
}

PerformExperiments <- function(parallel_ = TRUE, nin.chunks = 50){
  if (parallel_ == TRUE){
    registerDoMC(cores = 50)
  }
  # Setup H2O.
  h2o.init(nthreads = -1)
  h2o.no_progress()

  # Load gene list.
  all.genes <- LoadGeneList()

  genes.loadorder <- CreateGeneLoadOrder(all.genes, nin.chunks)
  metaplan.out <- CreateMetaGeneList(all.genes, nin.chunks)
  PerformRoundB(all.genes, metaplan.out, genes.loadorder) 
}

PerformExperiments()
