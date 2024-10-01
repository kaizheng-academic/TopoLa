# IRIS+TopoLa 
## overview of IRIS+TopoLa 

<p align="center">
<img src="https://github.com/kaizheng-academic/TopoLa/blob/main/src/IRIS_TopoLa.png" width="1000" />
</p>
Overview of the IRIS+TopoLa workflow. First, IRIS+TopoLa requires as input spatial transcriptomics (ST) data measured across multiple tissue slices with spatial localization information, along with scRNA-seq reference data from the same tissue containing cell-type-specific gene expression information. Second, IRIS+TopoLa integrates scRNA-seq data into the ST data and consolidates ST data across multiple tissue slices to perform comprehensive reference domain detection and obtain cell representations. Third, these cell representations are fed into the TopoLa framework to obtain enhanced cell representations. Finally, clustering is performed based on the enhanced cell representations to identify spatial domains within the tissue.   

Installation
------------

You need to install IRIS first. Click [here](https://github.com/YingMa0107/IRIS/)  for detailed installation instructions.


## Datasets 

You need to click [here](https://drive.google.com/file/d/1I4LIoFcd1kSP4SGpsNAMPumnROuz9N6A/view?usp=share_link)  to download the relevant dataset and place all the files into the `Data` folder.


How to use `demo.R`
-------------------
To directly use demo.R, you can run the following code:
```R
Rscript demo.R
```
