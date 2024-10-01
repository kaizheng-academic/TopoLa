## overview of scGPT+TopoLa for multi-batch integration

<p align="center">
<img src="https://github.com/kaizheng-academic/TopoLa/blob/main/src/scGPT_TopoLa_omic.png" width="1000" />
</p>
Enhance cell representations of scGPT for multi-omic integration using the TopoLa framework. The model begins with generative pretraining on large-scale scRNA-seq data from the cell atlases. Next, the pretrained model parameters are fine-tuned using multi-omics single-cell data. A cell graph, constructed from the learned cell embeddings, is then input into the TopoLa framework to obtain an enhanced cell graph. Finally, this enhanced cell graph is input into the Louvain algorithm for cell clustering.

Installation
------------

1. Install **scGPT** by following the instructions provided [here](https://github.com/bowang-lab/scGPT).
2. Download the **scGPT** code by clicking [here](https://github.com/bowang-lab/scGPT/archive/refs/heads/main.zip), and then unzip the files.
3. Copy the following files: `Louvain_TopoLa.py`, `demo.py`, `scGPT_TopoLa.py`, and `utils_TopoLa.py` into the `tutorials` folder.


## pretrained model 

To use the pretrained model, follow the steps below:

1.	Download the pretrained model by clicking [here](https://drive.google.com/file/d/1__se85Ru86rS7By4Zwvbel5VinEDPlmj/view?usp=share_link).
2.	In the `tutorials` folder, create a new folder named `save`.
3.	Inside the `save` folder, create another folder called `BMMC`.
4.	Unzip the downloaded BMMC.zip file.
5.	Move the extracted best_model.pt file into the `BMMC` folder.


How to use `demo.m`
-------------------
To directly use demo.m, you can run the following code:
```python
python demo.py
```
