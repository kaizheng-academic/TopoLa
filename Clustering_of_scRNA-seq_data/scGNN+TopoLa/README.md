# scGNN+TopoLa 
## overview of scGNN+TopoLa 

<p align="center">
<img src="https://github.com/kaizheng-academic/TopoLa/blob/main/src/scGNN_TopoLa.png" width="1000" />
</p>
Enhance cell representations of scGNN using the TopoLa framework. Firstly, the gene expression matrix is used as input. A feature autoencoder learns dimensional representations of the input as embeddings, which are subsequently used to construct and refine a cell graph. A graph autoencoder then learns the topological embeddings of the cell graph for clustering. Within each cell type, a distinct clustering autoencoder is used to reconstruct gene expression values. The framework iteratively uses the reconstructed expressions as new input until convergence is achieved. Subsequently, an imputation autoencoder regularizes the imputed gene expression values based on the cell-cell relationships learned from the cell graph, applying this regularization to the original preprocessed expression matrix. The final learned cell graph is used as the cell representation input for the TopoLa framework, where the enhanced cell graph is fed into the Louvain algorithm for cell clustering.

Installation
------------

Before running the “demo.m”, you need to install MATLAB and compile the two `.cpp` files. Please follow these steps:

1. Install MATLAB.
2. Compile the two `.cpp` files before running the Louvain method.

Run the following commands in MATLAB:

```matlab
mex jl_mnew.cpp
mex jl_clust.cpp
```

Once the compilation is complete, you can proceed with running the  “demo.m".


## Datasets 

You need to click [here](https://drive.google.com/file/d/1k3AD1tdELfYM1xXs9jjhKJ6nb1pu71As/view?usp=share_link)  to download the relevant dataset and click [here](https://drive.google.com/file/d/1JsN3IXCn8KzEZGuwMNALIBVdxSeybAI7/view?usp=share_link)  to download the result_scGNN.


How to use `demo.m`
-------------------
To directly use demo.m, you can run the following code:
```matlab
matlab -r demo
```
