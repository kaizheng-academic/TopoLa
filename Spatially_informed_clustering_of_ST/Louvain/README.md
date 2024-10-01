Installation
------------

Before running the “demo.m”, you need to install MATLAB and compile the two `.cpp` files. Please follow these steps:

1. Install MATLAB.
2. Compile the two `.cpp` files before running the Louvain method.

Run the following commands in MATLAB:

```matlab
mex jl_mnew.cpp
mex jl_clust.cpp
```

Once the compilation is complete, you can proceed with running the  “demo.m".


Datasets 
-------------------
You need to download "graphST_mculst_data.rar" from this [link](https://drive.google.com/file/d/1w1Ghtt7mq5qHvD6DQ-vJKLTKvZ-oxrpz/view?usp=sharing) and place all the files into the `data` folder.


How to use `demo.m`
-------------------
To directly use demo.m, you can run the following code:
```matlab
matlab -r demo
```
