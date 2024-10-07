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

You need to click [here](https://github.com/kaizheng-academic/TopoLa/tree/main/Rare_cell_identification/data)  to download the relevant dataset.


How to use `demo.py`
-------------------
To directly use demo.py, you can run the following code:
```python
python demo.py
```


## Dependencies

SCA requires the following packages:

absl-py                        1.2.0
anndata                        0.8.0
anyio                          3.6.1
argon2-cffi                    21.3.0
argon2-cffi-bindings           21.2.0
attrs                          22.1.0
Babel                          2.10.3
backcall                       0.2.0
beautifulsoup4                 4.11.1
bleach                         5.0.1
Bottleneck                     1.3.4
brotlipy                       0.7.0
cached-property                1.5.2
cachetools                     5.2.0
certifi                        2022.12.7
cffi                           1.14.6
chardet                        4.0.0
cmake                          3.24.1
conda                          23.1.0
conda-package-handling         1.7.3
cryptography                   3.4.7
cycler                         0.11.0
debugpy                        1.6.3
decorator                      5.1.1
defusedxml                     0.7.1
entrypoints                    0.4
fastjsonschema                 2.16.1
flit-core                      3.6.0
fonttools                      4.36.0
google-auth                    2.10.0
google-auth-oauthlib           0.4.6
grpcio                         1.47.0
h5py                           3.6.0
idna                           2.10
igraph                         0.10.8
importlib-metadata             4.12.0
importlib-resources            5.9.0
ipykernel                      6.15.1
ipython                        7.34.0
ipython-genutils               0.2.0
ipywidgets                     8.0.1
jedi                           0.18.1
Jinja2                         3.1.2
joblib                         1.3.2
json5                          0.9.10
jsonschema                     4.12.1
jupyter-client                 7.3.4
jupyter-core                   4.11.1
jupyter-server                 1.18.1
jupyterlab                     3.4.5
jupyterlab-language-pack-zh-CN 3.4.post3
jupyterlab-pygments            0.2.2
jupyterlab-server              2.15.0
jupyterlab-widgets             3.0.1
kiwisolver                     1.4.4
leidenalg                      0.10.2
llvmlite                       0.39.1
lxml                           4.9.1
Markdown                       3.4.1
MarkupSafe                     2.1.1
matplotlib                     3.5.3
matplotlib-inline              0.1.6
mistune                        0.8.4
mkl-fft                        1.3.1
mkl-random                     1.2.2
mkl-service                    2.4.0
natsort                        8.4.0
nbclassic                      0.4.3
nbclient                       0.6.6
nbconvert                      6.5.3
nbformat                       5.4.0
nest-asyncio                   1.5.5
networkx                       2.6.3
notebook                       6.4.12
notebook-shim                  0.1.0
numba                          0.56.4
numexpr                        2.7.3
numpy                          1.21.6
oauthlib                       3.2.0
packaging                      21.3
pandas                         1.3.5
pandocfilters                  1.5.0
parso                          0.8.3
patsy                          0.5.6
pexpect                        4.8.0
pickleshare                    0.7.5
Pillow                         9.2.0
pip                            21.1.3
pkgutil-resolve-name           1.3.10
pluggy                         1.0.0
prometheus-client              0.14.1
prompt-toolkit                 3.0.30
protobuf                       3.19.4
psutil                         5.9.1
ptyprocess                     0.7.0
pyasn1                         0.4.8
pyasn1-modules                 0.2.8
pycosat                        0.6.3
pycparser                      2.20
Pygments                       2.13.0
pynndescent                    0.5.13
pyOpenSSL                      20.0.1
pyparsing                      3.0.9
pyrsistent                     0.18.1
PySocks                        1.7.1
python-dateutil                2.8.2
python-igraph                  0.10.8
pytz                           2022.7
pyzmq                          23.2.1
requests                       2.25.1
requests-oauthlib              1.3.1
rsa                            4.9
ruamel-yaml-conda              0.15.100
ruamel.yaml                    0.16.12
ruamel.yaml.clib               0.2.6
scanpy                         1.9.3
scikit-learn                   1.0.2
scipy                          1.7.3
seaborn                        0.12.2
Send2Trash                     1.8.0
session-info                   1.0.0
setuptools                     52.0.0.post20210125
shannonca                      0.0.7
six                            1.16.0
sklearn                        0.0
sniffio                        1.2.0
soupsieve                      2.3.2.post1
statsmodels                    0.13.5
stdlib-list                    0.10.0
tensorboard                    2.10.0
tensorboard-data-server        0.6.1
tensorboard-plugin-wit         1.8.1
terminado                      0.15.0
texttable                      1.7.0
threadpoolctl                  3.1.0
tinycss2                       1.1.1
toolz                          0.12.0
tornado                        6.2
tqdm                           4.61.2
traitlets                      5.3.0
typing-extensions              4.4.0
umap-learn                     0.5.6
urllib3                        1.26.6
wcwidth                        0.2.5
webencodings                   0.5.1
websocket-client               1.3.3
Werkzeug                       2.2.2
wheel                          0.36.2
widgetsnbextension             4.0.1
zipp                           3.11.0
