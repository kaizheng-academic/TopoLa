# Deciphering Complex Networks via Topology-encoded Latent Hyperbolic Geometry
Complex networks, abstracting real-world systems, pose cross-disciplinary challenges in decoding their core information. Recently, latent hyperbolic geometry has gained traction in network analysis due to its ability to preserve nodes' local intrinsic properties. This study thoroughly reveals a deep and strong connection between global topological structures and latent space positioning, which prompts us to propose a novel embedding framework called Topology-encoded Latent Hyperbolic Geometry (TopoLa) for analyzing complex networks. With the encoded topological information in the latent space, TopoLa is capable of enhancing both conventional and low-rank networks, using the singular value gap to clarify the mathematical principles behind this enhancement. Meanwhile, we show that the equipped TopoLa distance can also help augment pivotal deep learning models encompassing knowledge distillation and contrastive learning.

# Requirements
* Matlab >= 2022

# Installation
LRTM can be downloaded by
```
git clone https://github.com/kaizheng-academic/TopoLa
```
Installation has been tested in a Windows platform.

# File Description
* link_prediction: "link prediction in complex networks" in the main test;
* singlecell: "Single-cell clustering via single-cell RNA-seq data" in the main test;
* HiC: "Domain identification via Hi-C networks;
* SpeciesIdentification: Fine-grained species identification;
* low-rank: Using NR for enhancing low-rank matrix;
* fastNR: disease ontology similarity matrix;
* didr: disease-drug association matrix.

# Functions Description
* ```LRTM.m```: this function can implement the LRTM algorithm;


# Instructions
We provide detailed step-by-step instructions for running LRTM model.

**Step 1**: add datasets\functions paths
```
addpath('Datasets');
addpath('Functions');
```
**Step 2**: load datasets with association matirx and similarity matrices
```
load Fdataset_ms
A_DR = didr;
R = (drug_AtcS+drug_TargetS)/2;
D = (disease_PhS+disease_DoS)/2;
```
**Step 3**: parameter Settings

The hyper-parameters are fixed.
```
alpha = 15; 
beta = 20; 

```

**Step 4**: run the LRTM algorithm
```
A_recovery = LRTM(A_DR',alpha,beta, R, D);
```

# A Quickstart Guide
Users can immediately start playing with LRTM running ```Demo_LRTM.m``` in matlab.
* ```Demo_LRTM.m```: it demonstrates a process of predicting drug-disease associations on Fdataset_ms by LRTM algorithm.

# Run LRTM on User's Own Data
We provided instructions on implementing LRTM model with user's own data. One could directly run LRTM model in ```Demo_LRTM.m``` with custom data by the following instructions.

**Step 1**: Prepare your own data and add into the ```Datasets``` folder

The required data includes drug-disease association matirx and similarity matrices, which are all saved by ```mat``` files.

**Step 2**: Modify four lines in ```Demo_LRTM.m```

You can find ```Fdataset_ms, A_DR, R, D``` in ```Demo_LRTM.m```. All you need to do is to replace them with your own data.



# Contact
If you have any questions or suggestions with the code, please let us know. 
Contact Kai Zheng at ```kaizheng@csu.edu.cn```
