# SCA+TopoLa 
## overview of SCA+TopoLa 

<p align="center">
<img src="https://github.com/kaizheng-academic/TopoLa/blob/main/src/SCA_TopoLa.png" width="1000" />
</p>
Enhance cell representations of SCA using the TopoLa framework. The model begins with the gene expression matrix as input. A surprisal score matrix is then generated using Smart-seq3. Singular value decomposition (SVD) is performed on the surprisal scores across all genes, producing right-eigenvectors that capture key axes of variation within the data. The input transcript count matrix is linearly projected onto these axes to generate a representation for each cell, serving as the initial cell representation. These cell representations are then input into the TopoLa framework to obtain enhanced representations, which are subsequently used for rare cell identification.    

Installation
------------

You need to install SCA first. Click [here](https://github.com/bendemeo/shannonca?tab=readme-ov-file)  for detailed installation instructions.


## Datasets 

You need to click [here](https://drive.google.com/file/d/1IEqq16XulvrirYL-oCGPbfRIzFz3OMS3/view?usp=share_link)  to download the relevant dataset.


How to use `demo.py`
-------------------
To directly use demo.py, you can run the following code:
```python
python demo.py
```


## Dependencies

SCA requires the following packages (python = 3.7.10):



scikit-learn = 1.0.2

scipy = 1.7.3

numpy = 1.21.6 

matplotlib

pandas 

seaborn

scanpy = 1.7.3

fbpca