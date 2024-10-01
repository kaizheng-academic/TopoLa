
## overview of graphST+TopoLa 

<p align="center">
<img src="https://github.com/kaizheng-academic/TopoLa/blob/main/src/graphST_TopoLa.png" width="1000" />
</p>
Overview of GraphST+TopoLa for spatial informed clustering. The model begins with preprocessed spatial gene expression data and a neighborhood graph constructed using spatial coordinates as inputs. A graph-based self-supervised contrastive learning approach is then employed to learn latent representations that preserve key features from gene expression profiles, spatial location, and local context information. These latent representations are mapped back to the original feature space to reconstruct the cell representation matrix. The matrix is subsequently input into the TopoLa framework to obtain refined cell representations. Finally, the cells are clustered based on their refined representations using the mclust algorithm.  

Installation
------------

You need to install graphST first. Click [here](https://github.com/bendemeo/shannonca?tab=readme-ov-file)  for detailed installation instructions.


## Datasets 

You need to click [here](https://github.com/kaizheng-academic/TopoLa/tree/main/Rare_cell_identification/data)  to download the relevant dataset.


How to use `demo.m`
-------------------
To directly use TopoLa.m, you can run the following code:
```python
python demo.py
```
