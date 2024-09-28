import pandas as pd
import scanpy as sc
import numpy as np
import os
import h5py
import warnings
from sklearn.metrics import f1_score
from collections import Counter
import time
import src.Leiden_TopoLa

warnings.filterwarnings("ignore")

def normalize(adata):
    sc.pp.filter_genes(adata, min_cells=3)
    adata.raw = adata.copy()
    sc.pp.normalize_per_cell(adata)
    adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
    sc.pp.log1p(adata)
    return adata

def calculate_rare_cell_types(y, res, h, n_cells):
    rare_types = [label for label, count in Counter(y).items() if count / n_cells < h]
    y_rare = [1 if label in rare_types else 0 for label in y]
    
    rare_clusters = [label for label, count in Counter(res).items() if count / n_cells < h]
    y_pred = [1 if label in rare_clusters else 0 for label in res]
    
    return y_rare, y_pred

def process_dataset(dataName, dir='./data'):
    print(f"Processing {dataName}")
    
    if dataName in ['10X_PBMC', 'Cao', 'Macosko', 'MacParland']:
        data_mat = h5py.File(os.path.join(dir, dataName + '.h5'))
        X = np.array(data_mat['X'])
        y = np.array(data_mat['Y'])
        geneName = ['gene'+str(i) for i in range(X.shape[1])]
        cellName = ['cell'+str(i) for i in range(X.shape[0])]
        data_mat.close()
        
    elif dataName in ["Chung", "Darmanis", "Deng", "Goolam", "Koh", "Li", "Pollen", "Puram", "Yang", "Zelsel"]:
        counts_df = pd.read_csv(os.path.join(dir, dataName + '.csv'), encoding='gbk')
        lab_df = pd.read_csv(os.path.join(dir, dataName + '_label.csv'))
        geneName = counts_df.iloc[:, 0].values.tolist()
        cellName = counts_df.columns[1:].values.tolist()
        y = lab_df.iloc[:, 1].values.tolist()
        X = counts_df.iloc[:, 1:].values.T
        y = np.array(y)
        
    elif dataName in ['Airway', 'hippocampus', 'iLNs', 'livers', 'Tonsil', 'UUOkidney',  'pancrea']:
        data_mat = h5py.File(os.path.join(dir, dataName + '.h5'))
        X = np.array(data_mat['X'])
        y = np.array(data_mat['Y'])
        geneName = np.array(data_mat['gn'])
        cellName = np.array(data_mat['cn'])
        data_mat.close()
        y = np.array([str(i, 'UTF-8') for i in y])
        geneName = np.array([str(i, 'UTF-8') for i in geneName])
        cellName = np.array([str(i, 'UTF-8') for i in cellName])
        


    n_cells = len(cellName)
    adata = sc.AnnData(X)
    adata.var_names = geneName
    adata.obs_names = cellName
    adata = normalize(adata)

    # Dimensionality reduction (using SCA)
    from shannonca.dimred import reduce_scanpy
    start_time = time.perf_counter()
    reduce_scanpy(adata, n_comps=50, iters=5, keep_scores=True, keep_loadings=True)
    end_time = time.perf_counter()
    print("Reduction time: {:.3f} seconds".format(end_time - start_time))

    # Neighbors calculation with error handling
    try:
        sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_sca', n_pcs=50)
    except Exception as e:
        print(f"Error in sc.pp.neighbors: {e}. Retrying...")
        sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_sca', n_pcs=50)

    # Fixed parameters
    h = 0.05
    alpha = 1e-3

    sc.tl.leiden(adata, TopoLa_bool=True, TopoLa_alpha=alpha)
    TopoLa_res = list(adata.obs['leiden'])
    
    y_rare, y_pred_TopoLa = calculate_rare_cell_types(y, TopoLa_res, h, n_cells)
    score2 = f1_score(y_rare, y_pred_TopoLa)
    print(f"F1_score for alpha {alpha} and h {h}: {score2}")
    
    return {
        "dataset": dataName,
        "alpha": alpha,
        "h": h,
        "F1_score": score2
    }

if __name__ == "__main__":
    datasets = ['Airway', 'UUOkidney', 'MacParland', '10X_PBMC', 'Cao', 'Chung', 
                'Darmanis', 'Deng', 'Goolam', 'Koh', 'Li', 'Pollen', 
                'Puram', 'Yang', 'Zelsel', 'hippocampus', 'iLNs', 
                'livers', 'Tonsil', 'pancrea']
    
    all_results = []
    
    for dataName in datasets:
        result = process_dataset(dataName)
        all_results.append(result)
    
    results_df = pd.DataFrame(all_results)
    results_df.to_csv("Result.csv", index=False)
    print("All results saved to Result.csv.")
