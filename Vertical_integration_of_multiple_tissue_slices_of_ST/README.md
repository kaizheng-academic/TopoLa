
## overview of graphST+TopoLa for multi-batch integration

<p align="center">
<img src="https://github.com/kaizheng-academic/TopoLa/blob/main/src/graphST_TopoLa_batch.png" width="1000" />
</p>
Overview of GraphST+TopoLa for multi-batch integration. First, hematoxylin and eosin (H&E) stained images from two or more samples are aligned to establish spatial correspondence across samples. A shared neighborhood graph is constructed, accounting for both intra- and inter-sample neighbors, which enables feature smoothing across the entire dataset. Second, GraphST+TopoLa corrects batch effects by smoothing features across samples to generate cell representations. Third, these cell representations are input into the TopoLa framework to obtain enhanced representations for further batch effect correction.

Installation
------------

You need to install graphST first. Click [here](https://github.com/JinmiaoChenLab/GraphST)  for detailed installation instructions.


## Datasets 

You need to click [here](https://drive.google.com/file/d/1IEqq16XulvrirYL-oCGPbfRIzFz3OMS3/view?usp=share_link)  to download the relevant dataset and place “filtered_feature_bc_matrix.h5ad” and “metadata_sample1_section2.tsv” into the `data` folder. Additionally, place “GraphST_VerticalST_final.h5ad” in the root directory.



