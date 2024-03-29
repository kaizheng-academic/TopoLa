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
* degree_TopoSimilarity: "Comparing the two energy distance measures" in the main test.
* link_prediction: "link prediction in complex networks" in the main test;
* singlecell: "Single-cell clustering via single-cell RNA-seq data" in the main test;
* HiC: "Domain identification via Hi-C networks“ in the main test;
* SpeciesIdentification: "Fine-grained species identification" in the main test;
* low-rank: "Using NR for enhancing low-rank matrix" in the main test;
* RKD: "Knowledge distillation" in the main test;
* MoCo: "Contrastive learning" in the main test;
* fastNR: "Comparison of performance between NR and fastNR" in the supplementary materials;



# Instructions
We provide the experimental code and execution steps involved in each section.

**Comparing the two energy distance measures**: The properties of energy distance measures.
```
matlab -r Neighbor_Dtopola
matlab -r TopoSimilarity_Dtopola
```
**link prediction in complex networks**: The performance of TRWR in link prediction across nine distinct complex networks.
```
matlab -r demo
```
**Single-cell clustering via single-cell RNA-seq data**: Identification of cell types using single-cell RNA-seq data

```
matlab -r demo 

```

**Domain identification via Hi-C networks**: Performance evaluation of boundary detection
```
matlab -r demo


```
**Fine-grained species identification**: Average species identification accuracy.
```
matlab -r demo 
```

**Using NR for enhancing low-rank matrix**: Utilizing low-rank matrix completion for movie recommendations and multi-label learning. 
```
matlab -r demo_Movielens

matlab -r demo_Multi_Label 
```

**Knowledge distillation**: The performance of RKD-Topola, their teacher (baseline) and student (baseline) models. 
```
python run.py --help    
python run_distill.py --help

# Train a teacher embedding network of resnet50 (d=512)
# using triplet loss (margin=0.2) with distance weighted sampling.
python run.py --mode train \ 
               --dataset cub200 \
               --base resnet50 \
               --sample distance \ 
               --margin 0.2 \ 
               --embedding_size 512 \
               --save_dir teacher

# Evaluate the teacher embedding network
python run.py --mode eval \ 
               --dataset cub200 \
               --base resnet50 \
               --embedding_size 512 \
               --load teacher/best.pth 

# Distill the teacher to student embedding network
python run_distill.py --dataset cub200 \
                      --base resnet18 \
                      --embedding_size 64 \
                      --l2normalize false \
                      --teacher_base resnet50 \
                      --teacher_embedding_size 512 \
                      --teacher_load teacher/best.pth \
                      --dist_ratio 1  \
                      --angle_ratio 2 \
                      --save_dir student
                      
# Distill the trained model to student network
python run.py --mode eval \ 
               --dataset cub200 \
               --base resnet18 \
               --l2normalize false \
               --embedding_size 64 \
               --load student/best.pth 
            
```

# Contact
If you have any questions or suggestions with the code, please let us know. 
Contact Kai Zheng at ```kaizheng@csu.edu.cn```
