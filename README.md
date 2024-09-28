# TopoLa
## TopoLa: A Universal Framework to Enhance Cell Representations for Single-cell and Spatial Omics through Topology-encoded Latent Hyperbolic Geometry

<p align="center">
<img src="https://github.com/kaizheng-academic/TopoLa/blob/main/src/TopoLa.png" width="400" />
</p>
Recent advances in cellular research demonstrate that scRNA-seq characterizes cellular heterogeneity, while spatial transcriptomics reveals the spatial distribution of gene expression. Cell representation is the fundamental issue in the two fields. Here, we propose Topology-encoded Latent Hyperbolic Geometry (TopoLa), a computational framework enhancing cell representations by capturing fine-grained intercellular topological relationships. The framework introduces a new metric, TopoLa distance (TLd), which quantifies the geometric distance between cells within latent hyperbolic space, capturing the network’s topological structure more effectively. With this framework, the cell representation can be enhanced considerably by performing convolution on its neighboring cells. Performance evaluation across seven biological tasks, including scRNA-seq data clustering and spatial transcriptomics domain identification, shows that TopoLa significantly improves the performance of several state-of-the-art models. These results underscore the generalizability and robustness of TopoLa, establishing it as a valuable tool for advancing both biological discovery and computational methodologies.


Installation
------------
You can install the released version of IRIS from Github with the following code, for more installation details or solutions that might solve related issues (specifically MacOS system) see the [link](https://yingma0107.github.io/IRIS/documentation/02_installation.html).

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
