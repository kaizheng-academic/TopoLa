# scGNN+TopoLa 
## overview of scGNN+TopoLa 

<p align="center">
<img src="https://github.com/kaizheng-academic/TopoLa/blob/main/src/SIMLR_TopoLa.png" width="1000" />
</p>
Enhance cell representations of SIMLR using the TopoLa framework. The pipeline begins with a gene expression matrix as input. SIMLR+TopoLa calculates multiple kernels to learn intercellular similarities. This similarity matrix is then input into the TopoLa framework, refining cellular similarities. These similarities are subsequently processed through dimensionality reduction techniques to generate cellular representations. Finally, cell clustering is performed using the affinity propagation algorithm.

Installation
------------

Before running the “demo.m”, you need to install MATLAB.


## Datasets 

You need to click [here](https://github.com/kaizheng-academic/TopoLa/tree/main/Clustering_of_scRNA-seq_data/scGNN%2BTopoLa/data)  to download the relevant dataset.


How to use `demo.m`
-------------------
To directly use demo.m, you can run the following code:
```matlab
matlab -r demo
```
