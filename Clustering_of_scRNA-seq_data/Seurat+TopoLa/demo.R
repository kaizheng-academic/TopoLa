suppressPackageStartupMessages({
  library(R.matlab)
  library(Seurat)
  library(Matrix)
  library(aricode)
})

set.seed(42)

target_file <- "./data/Baron_human1.mat" 
fixed_lambda <- 1e-6

TopoLa_Graph_Enhance <- function(snn_matrix, lambda) {
  A <- as.matrix(snn_matrix)
  I <- diag(nrow(A))
  AtA <- t(A) %*% A
  
  C <- tryCatch({
    solve((1/lambda) * I + AtA, AtA)
  }, error = function(e) {
    solve((1/lambda) * I + AtA , AtA)
  })
  
  A_enhanced <- A %*% C
  
  A_enhanced[is.na(A_enhanced)] <- 0
  A_enhanced[is.infinite(A_enhanced)] <- 0
  A_enhanced[A_enhanced < 0] <- 0
  diag(A_enhanced) <- 0
  
  rownames(A_enhanced) <- rownames(snn_matrix)
  colnames(A_enhanced) <- colnames(snn_matrix)
  return(A_enhanced)
}

m <- readMat(target_file)
counts <- if (!is.null(m$in.X)) m$in.X else m$in_X
counts <- as(counts, "dgCMatrix") 

raw_labs <- as.numeric(m$labs)
type_names <- as.character(unlist(m$id))
true_types <- if (min(raw_labs) == 0) type_names[raw_labs + 1] else type_names[raw_labs]
true_num <- as.numeric(as.factor(true_types))

if (ncol(counts) != length(raw_labs)) counts <- t(counts)
colnames(counts) <- paste0("C", 1:ncol(counts))
rownames(counts) <- paste0("G", 1:nrow(counts))

so <- CreateSeuratObject(counts = counts)
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)

safe_npcs <- min(30, (ncol(so)-1), (nrow(so)-1))
so <- RunPCA(so, npcs = safe_npcs, verbose = FALSE, seed.use = 42)
so <- FindNeighbors(so, dims = 1:min(20, safe_npcs), verbose = FALSE)

enhanced_mat <- TopoLa_Graph_Enhance(so[["RNA_snn"]], fixed_lambda)
enhanced_graph <- as.Graph(as(enhanced_mat, "dgCMatrix"))
DefaultAssay(enhanced_graph) <- "RNA"
so[["topo_snn"]] <- enhanced_graph

so <- FindClusters(so, graph.name = "topo_snn", resolution = 0.5, verbose = FALSE, random.seed = 42)

ari_score <- ARI(as.numeric(so$seurat_clusters), true_num)

cat(sprintf("Dataset: %s | Lambda: %s | ARI: %.4f\n", basename(target_file), fixed_lambda, ari_score))