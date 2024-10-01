
source("function.R")

# Loop through D values from 1 to 3
for (D in 1:3) {
  
  # Load corresponding data based on D value
  if (D == 1) {
    load(file = "./Data/07-10_MatrixCombined.Rdata")
  } else if (D == 2) {
    load(file = "./Data/69-72_MatrixCombined.Rdata")
  } else if (D == 3) {
    load(file = "./Data/73-76_MatrixCombined.Rdata")
  }
  
  # Parameters
  lambda_values <- c(1e4)
  datasets <- unique(MatrixCombined$Slice)
  
  results <- data.frame(Dataset = character(), lambda = numeric(), NMI = numeric(), RI = numeric(), ARI = numeric(), stringsAsFactors = FALSE)
  
  for (i in datasets) {
    best_ari <- -Inf
    best_lambda <- NULL
    best_nmi <- NULL
    best_ri <- NULL
    best_clusters <- NULL
    
    slice_data <- MatrixCombined[MatrixCombined$Slice == i, ]
    
    for (lambda in lambda_values) {
      # Apply TopoLa function
      MatrixProcessed <- TopoLa(slice_data[, CommonCellTypes], lambda)
      
 
      
      # K-means clustering
      kmeans <- kmeansFunc_Iter(MatrixProcessed, 5, centers_old)
      numCluster <- length(unique(kmeans$kmeans))
      slice_data$kmeans_cluster <- as.numeric(kmeans$kmeans)
      slice_data$kmeans_cluster <- slice_data$kmeans_cluster - 1
      
      # Evaluate clustering
      lab <- read.table(paste0("./Data/", i, "_labs.txt"))[, 1]
      iris <- slice_data$kmeans_cluster
      metrics <- evalcluster(lab, iris)
      print(paste("lambda:", lambda, "Slice:", i, "Results:", metrics["ARI"]))
      
      if (metrics["ARI"] > best_ari) {
        best_ari <- metrics["ARI"]
        best_lambda <- lambda
        best_nmi <- metrics["NMI"]
        best_ri <- metrics["RI"]
        best_clusters <- slice_data$kmeans_cluster
      }
    }
    
    # Save results to the results data frame
    results <- rbind(results, data.frame(Dataset = i, lambda = best_lambda, NMI = best_nmi,  ARI = best_ari))
  }
  
  # Print results
  print(results)
  
  # Append results to Result.csv
  result_file <- "./Result.csv"
  
  # Check if file exists to determine whether to write header
  if (file.exists(result_file)) {
    write.table(results, file = result_file, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  } else {
    write.csv(results, file = result_file, row.names = FALSE)
  }
}