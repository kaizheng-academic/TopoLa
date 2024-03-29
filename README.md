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
**Single-cell clustering via single-cell RNA-seq data**: Identification of cell types using single-cell RNA-seq data (Please decompress the file before executing)

```
matlab -r demo 
```

**Domain identification via Hi-C networks**: Performance evaluation of boundary detection (Please decompress the file before executing)
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
# Python=3.6,Pytorch=1.0,tqdm,h5py,scipy 
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

**Contrastive learning**: The performance of MoCo-Topola. 
Usage: Preparation

Install PyTorch and download the ImageNet dataset following the [official PyTorch ImageNet training code](https://github.com/pytorch/examples/tree/master/imagenet). Similar to [MoCo v1/2](https://github.com/facebookresearch/moco), this repo contains minimal modifications on the official PyTorch ImageNet code. We assume the user can successfully run the official PyTorch ImageNet code.
For ViT models, install [timm](https://github.com/rwightman/pytorch-image-models) (`timm==0.4.9`).

The code has been tested with CUDA 10.2/CuDNN 7.6.5, PyTorch 1.9.0 and timm 0.4.9.

Usage: Self-supervised Pre-Training

Below are three examples for MoCo v3 pre-training. 

ResNet-50 with 2-node (16-GPU) training, batch 4096

On the first node, run:
```
python main_moco.py \
  --moco-m-cos --crop-min=.2 \
  --dist-url 'tcp://[your first node address]:[specified port]' \
  --multiprocessing-distributed --world-size 2 --rank 0 \
  [your imagenet-folder with train and val folders]
```
On the second node, run the same command with `--rank 1`.
With a batch size of 4096, the training can fit into 2 nodes with a total of 16 Volta 32G GPUs. 


ViT-Small with 1-node (8-GPU) training, batch 1024

```
python main_moco.py \
  -a vit_small -b 1024 \
  --optimizer=adamw --lr=1.5e-4 --weight-decay=.1 \
  --epochs=300 --warmup-epochs=40 \
  --stop-grad-conv1 --moco-m-cos --moco-t=.2 \
  --dist-url 'tcp://localhost:10001' \
  --multiprocessing-distributed --world-size 1 --rank 0 \
  [your imagenet-folder with train and val folders]
```

ViT-Base with 8-node training, batch 4096

With a batch size of 4096, ViT-Base is trained with 8 nodes:
```
python main_moco.py \
  -a vit_base \
  --optimizer=adamw --lr=1.5e-4 --weight-decay=.1 \
  --epochs=300 --warmup-epochs=40 \
  --stop-grad-conv1 --moco-m-cos --moco-t=.2 \
  --dist-url 'tcp://[your first node address]:[specified port]' \
  --multiprocessing-distributed --world-size 8 --rank 0 \
  [your imagenet-folder with train and val folders]
```
On other nodes, run the same command with `--rank 1`, ..., `--rank 7` respectively.

Notes:
1. The batch size specified by `-b` is the total batch size across all GPUs.
1. The learning rate specified by `--lr` is the *base* lr, and is adjusted by the [linear lr scaling rule](https://arxiv.org/abs/1706.02677) in [this line](https://github.com/facebookresearch/moco-v3/blob/main/main_moco.py#L213).
1. Using a smaller batch size has a more stable result (see paper), but has lower speed. Using a large batch size is critical for good speed in TPUs (as we did in the paper).
1. In this repo, only *multi-gpu*, *DistributedDataParallel* training is supported; single-gpu or DataParallel training is not supported. This code is improved to better suit the *multi-node* setting, and by default uses automatic *mixed-precision* for pre-training.

Usage: Linear Classification

By default, we use momentum-SGD and a batch size of 1024 for linear classification on frozen features/weights. This can be done with a single 8-GPU node.

```
python main_lincls.py \
  -a [architecture] --lr [learning rate] \
  --dist-url 'tcp://localhost:10001' \
  --multiprocessing-distributed --world-size 1 --rank 0 \
  --pretrained [your checkpoint path]/[your checkpoint file].pth.tar \
  [your imagenet-folder with train and val folders]
```

Usage: End-to-End Fine-tuning ViT

To perform end-to-end fine-tuning for ViT, use our script to convert the pre-trained ViT checkpoint to [DEiT](https://github.com/facebookresearch/deit) format:
```
python convert_to_deit.py \
  --input [your checkpoint path]/[your checkpoint file].pth.tar \
  --output [target checkpoint file].pth
```
Then run the training (in the DeiT repo) with the converted checkpoint:
```
python $DEIT_DIR/main.py \
  --resume [target checkpoint file].pth \
  --epochs 150
```

Usage: KNN
```
python main_knn.py
```

Note:
1. We use `--resume` rather than `--finetune` in the DeiT repo, as its `--finetune` option trains under eval mode. When loading the pre-trained model, revise `model_without_ddp.load_state_dict(checkpoint['model'])` with `strict=False`.
1. Our ViT-Small is with `heads=12` in the Transformer block, while by default in DeiT it is `heads=6`. Please modify the DeiT code accordingly when fine-tuning our ViT-Small model. 


# Contact
If you have any questions or suggestions with the code, please let us know. 
Contact Kai Zheng at ```kaizheng@csu.edu.cn```
