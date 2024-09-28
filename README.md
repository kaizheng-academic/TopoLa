# TopoLa
## TopoLa: A Universal Framework to Enhance Cell Representations for Single-cell and Spatial Omics through Topology-encoded Latent Hyperbolic Geometry

<p align="center">
<img src="https://github.com/kaizheng-academic/TopoLa/blob/main/src/TopoLa.png" width="400" />
</p>
Recent advances in cellular research demonstrate that scRNA-seq characterizes cellular heterogeneity, while spatial transcriptomics reveals the spatial distribution of gene expression. Cell representation is the fundamental issue in the two fields. Here, we propose Topology-encoded Latent Hyperbolic Geometry (TopoLa), a computational framework enhancing cell representations by capturing fine-grained intercellular topological relationships. The framework introduces a new metric, TopoLa distance (TLd), which quantifies the geometric distance between cells within latent hyperbolic space, capturing the network’s topological structure more effectively. With this framework, the cell representation can be enhanced considerably by performing convolution on its neighboring cells. Performance evaluation across seven biological tasks, including scRNA-seq data clustering and spatial transcriptomics domain identification, shows that TopoLa significantly improves the performance of several state-of-the-art models. These results underscore the generalizability and robustness of TopoLa, establishing it as a valuable tool for advancing both biological discovery and computational methodologies.




## Description 
* The code "TopoLa.m" is part of the TopoLa framework, implemented using MATLAB. To run this code, you need to have MATLAB installed, with no specific version requirement.
* The “Clustering_of_scRNA-seq_data” folder contains resources related to the task of clustering scRNA-seq data, including data, code, and tutorials.
* The "Rare_cell_identification” folder contains resources related to the task of rare cell identification, including data, code, and tutorials.
* The “Single_cell_multi_batch_integration” folder contains resources related to the task of single cell multi-batch integration, including data, code, and tutorials.
* The “Single_cell_multi_omics_integration” folder contains resources related to the task of single cell multi-omics integration, including data, code, and tutorials.
* The “Spatially_informed_clustering_of_ST” folder contains resources related to the task of spatially informed clustering of ST, including data, code, and tutorials.
* The “Cspatially_informed_clustering_with_the_integration_of_scRNA-seq_and_ST” folder contains resources related to the task of spatially informed clustering with the integration of scRNA-seq and ST, ncluding data, code, code, and tutorials.
* The “Vertical_integration_of_multiple_tissue_slices_of_ST” folder contains resources related to the task of spatially informed clustering with the integration of scRNA-seq and ST, including data, code, and tutorials.


The "TopoLa.m" has been run successfully on Operating systems: 
* macOS Sonoma 14.1.2
* Ubuntu 18.04.6 LTS
* Windows 10

## How to use `TopoLa`
-------------------
To directly use TopoLa.m, you can run the following code:
```matlab
matlab -r TopoLa
```
You need to provide the “input.mat” file in the directory where “TopoLa.m” is being run to store your matrix. After execution, the enhanced network will be exported as “output.mat” in the current directory. Click [here](https://drive.google.com/file/d/1Cl9BmLQS7zJ8SlkF1OxHuargU5fbUNFi/view?usp=share_link) to download an example of “input.mat”.


Examples are provided for seven key biological tasks, with details as follows:
| Biological research       | Related tutorials                                        | Dataset                                                                                     |
| :------------------------ | :------------------------------------------------------ | :------------------------------------------------------------------------------------------- |
| Single-cell multi-batch integration      | [scGPT](https://github.com/kaizheng-academic/TopoLa/tree/main/Single_cell_multi_batch_integration)             | [link](https://drive.google.com/drive/folders/1_GROJTzXiAV8HB4imruOTk6PEGuNOcgB?usp=sharing) |
| clustering scRNA-seq data | [link](https://github.com/kaizheng-academic/TopoLa/tree/main/Clustering_of_scRNA-seq_data/)   | [link](https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y?usp=sharing) |
| Single-cell multi-batch integration      | [link](https://github.com/kaizheng-academic/TopoLa/tree/main/Single_cell_multi_batch_integration)             | [link](https://drive.google.com/drive/folders/1_GROJTzXiAV8HB4imruOTk6PEGuNOcgB?usp=sharing) |
| Single-cell multi-omic integration    | [link](https://github.com/kaizheng-academic/TopoLa/tree/main/Single_cell_multi_omics_integration)                 | [link](https://drive.google.com/drive/folders/1vf1ijfQSk7rGdDGpBntR5bi5g6gNt-Gx?usp=sharing) |
| Rare cell identification      | [link](https://github.com/kaizheng-academic/TopoLa/tree/main/Rare_cell_identification)  | [link](https://drive.google.com/drive/folders/1kkug5C7NjvXIwQGGaGoqXTk_Lb_pDrBU?usp=sharing) |
| Spatially informed clustering of ST    | [link](https://github.com/kaizheng-academic/TopoLa/tree/main/Spatially_informed_clustering_of_ST)         | [link](https://drive.google.com/drive/folders/1GcgXrd7apn6y4Ze_iSCncskX3UsWPY2r?usp=sharing) |
| Vertical integration of multiple tissue slices  | [link](https://github.com/kaizheng-academic/TopoLa/tree/main/Vertical_integration_of_multiple_tissue_slices_of_ST)      | [link](https://drive.google.com/drive/folders/16A1DJ30PT6bodt4bWLa4hpS7gbWZQFBG?usp=sharing) |
| Spatially informed clustering with the integration of scRNA-seq and ST | [IRIS](https://github.com/kaizheng-academic/TopoLa/tree/main/spatially_informed_clustering_with_the_integration_of_scRNA-seq_and_ST)  | [link](https://drive.google.com/drive/folders/1S-1AR65DF120kNFpEbWCvRHPhpkGK3kK?usp=sharing) |



# Issues
We warmly welcome all feedback, bug reports, and suggestions! When raising issues, please ensure they include a detailed and reproducible example for easier troubleshooting.






