
## overview of graphST+TopoLa 

<p align="center">
<img src="https://github.com/kaizheng-academic/TopoLa/blob/main/src/graphST_TopoLa.png" width="1000" />
</p>
Overview of GraphST+TopoLa for spatial informed clustering. The model begins with preprocessed spatial gene expression data and a neighborhood graph constructed using spatial coordinates as inputs. A graph-based self-supervised contrastive learning approach is then employed to learn latent representations that preserve key features from gene expression profiles, spatial location, and local context information. These latent representations are mapped back to the original feature space to reconstruct the cell representation matrix. The matrix is subsequently input into the TopoLa framework to obtain refined cell representations. Finally, the cells are clustered based on their refined representations using the mclust algorithm.  

Installation
------------

You need to install graphST first. Click [here](https://drive.google.com/file/d/1w1Ghtt7mq5qHvD6DQ-vJKLTKvZ-oxrpz/view?usp=sharing)  for detailed installation instructions.


## Datasets 

You need to click [here](https://drive.google.com/file/d/1w1Ghtt7mq5qHvD6DQ-vJKLTKvZ-oxrpz/view?usp=sharing)  to download the relevant dataset.


How to run `graphST+TopoLa`
-------------------
In our experiments, we evaluated the performance of scGPT+TopoLa with two clustering methods:

1.	To view the tutorial for the Louvain algorithm, click [here](https://github.com/kaizheng-academic/TopoLa/tree/main/Spatially_informed_clustering_of_ST/Louvain).
2.	To view the tutorial for the mclust algorithm, click [here](https://github.com/kaizheng-academic/TopoLa/tree/main/Spatially_informed_clustering_of_ST/mclust).
