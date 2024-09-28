# scGNN+TopoLa 
## overview of scGNN+TopoLa 

<p align="center">
<img src="https://github.com/kaizheng-academic/TopoLa/blob/main/src/scGNN_TopoLa.png" width="1000" />
</p>
Enhance cell representations of scGNN using the TopoLa framework. Firstly, the gene expression matrix is used as input. A feature autoencoder learns dimensional representations of the input as embeddings, which are subsequently used to construct and refine a cell graph. A graph autoencoder then learns the topological embeddings of the cell graph for clustering. Within each cell type, a distinct clustering autoencoder is used to reconstruct gene expression values. The framework iteratively uses the reconstructed expressions as new input until convergence is achieved. Subsequently, an imputation autoencoder regularizes the imputed gene expression values based on the cell-cell relationships learned from the cell graph, applying this regularization to the original preprocessed expression matrix. The final learned cell graph is used as the cell representation input for the TopoLa framework, where the enhanced cell graph is fed into the Louvain algorithm for cell clustering.

Installation
------------

Before running the “demmo.m”, you need to install MATLAB and compile the two `.cpp` files. Please follow these steps:

1. Install MATLAB.
2. Compile the two `.cpp` files before running the Louvain method.

Run the following commands in MATLAB:

```matlab
mex jl_mnew.cpp
mex jl_clust.cpp
```

Once the compilation is complete, you can proceed with running the  “demmo.m.


## Dependencies 
* R version >= 4.2.2.
* Dependent R packages: Rcpp (>= 1.0.9), RcppArmadillo, SingleCellExperiment, SummarizedExperiment, methods, Matrix, MCMCpack, fields, wrMisc, RANN, stats, ggplot2, grDevices, reshape2


``` r
# install devtools if necessary
install.packages('devtools')

# install the IRIS package
devtools::install_github('YingMa0107/IRIS')

# load package
library(IRIS)

```
The R package has been installed successfully on Operating systems: 
* macOS Catalina 10.15, macOS Monterey 12.4
* Ubuntu 18.04.6 LTS
* Windows 10

# Issues
All feedback, bug reports and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible example and also please provide the output of your sessionInfo() in R! 

How to cite `IRIS`
-------------------
Ying Ma and Xiang Zhou. Integrative and Reference-Informed Spatial Domain Detection for Spatial Transcriptomics, Nature Methods 2024， [link](https://www.nature.com/articles/s41592-024-02284-9)

How to use `IRIS`
-------------------
Details in [Tutorial](https://yingma0107.github.io/IRIS/)
