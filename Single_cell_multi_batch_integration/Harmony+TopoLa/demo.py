import os

# 1. Threading Setup
THREADS = "1"
os.environ.update({
    "OMP_NUM_THREADS": THREADS, "MKL_NUM_THREADS": THREADS, 
    "OPENBLAS_NUM_THREADS": THREADS, "NUMEXPR_NUM_THREADS": THREADS,
    "OMP_PROC_BIND": "true", "OMP_PLACES": "cores"
})
print(f"Threads set to {THREADS}. Please Restart Kernel now.")
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import harmonypy as hm
import scib
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import svds

# --- 1. Global Settings ---
THREADS = "1"  # Set to your core count
os.environ.update({
    "OMP_NUM_THREADS": THREADS, "MKL_NUM_THREADS": THREADS, 
    "OPENBLAS_NUM_THREADS": THREADS, "NUMEXPR_NUM_THREADS": THREADS
})

SEED = 42
np.random.seed(SEED)
sc.set_figure_params(figsize=(4, 4), frameon=False)
sc.settings.verbosity = 1
try: sc.settings.n_jobs = int(THREADS)
except: pass

# --- 2. Data Loading & Preprocessing ---
print(">>> Loading and Preprocessing...")
adata = scvi.data.pbmc_dataset()
adata.obs["celltype"] = adata.obs["str_labels"].astype("category")
adata.obs["str_batch"] = adata.obs["batch"].astype(str).astype("category")
adata.var.set_index("gene_symbols", inplace=True)

# Fair Preprocessing (HVG not batch-aware)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1200, flavor="seurat_v3", subset=True)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30, svd_solver="randomized", random_state=SEED)

# Harmony Integration (Safe Shape Fix)
ho = hm.run_harmony(
    np.array(adata.obsm["X_pca"]), adata.obs, ["str_batch"], 
    max_iter_harmony=20, random_state=SEED
)
Z_harmony = np.array(ho.Z_corr)
if Z_harmony.shape[0] != adata.n_obs: 
    Z_harmony = Z_harmony.T
adata.obsm["X_pca_harmony"] = Z_harmony

# Neighbors & UMAP
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, random_state=SEED)
sc.tl.umap(adata, min_dist=0.3, random_state=SEED)

# --- 3. TopoLa Logic (Single Run) ---
class TopoLaHandler:
    def __init__(self, adata, k=1000):
        self.adata = adata
        self.adj = adata.obsp["connectivities"].tocsr()
        self.cache = None
        self.k = k

    def run_louvain(self, res, key_added, alpha=None, eps=-1e-6):
        adj_use = self.adj
        # Apply TopoLa filter if alpha (lambda) is provided
        if alpha is not None:
            if self.cache is None:
                A = (self.adj + self.adj.T) * 0.5
                A.eliminate_zeros()
                k_use = min(self.k, A.shape[0] - 2)
                U, s, VT = svds(A, k=k_use, which="LM")
                idx = np.argsort(s)[::-1]
                rows, cols = A.nonzero()
                self.cache = (U[:, idx], s[idx], VT[idx, :].T, rows, cols, A.shape)
            
            U, s, VT_T, rows, cols, shape = self.cache
            # TopoLa formula: s_new = s^3 / (s^2 + 1/alpha)
            s_new = (s**3) / (s**2 + (1.0 / alpha))
            
            # Reconstruct sparse matrix
            vals = np.sum((U[rows, :] * s_new) * VT_T[cols, :], axis=1)
            vals = np.maximum(vals, 0.0)
            
            vmin, vmax = vals.min(), vals.max()
            if vmax > vmin: vals = (vals - vmin) / (vmax - vmin)
            else: vals = np.zeros_like(vals)
            
            mask = vals >= eps
            W = csr_matrix((vals[mask], (rows[mask], cols[mask])), shape=shape)
            adj_use = (W + W.T) * 0.5
            adj_use.eliminate_zeros()

        sc.tl.louvain(
            self.adata, resolution=res, key_added=key_added, 
            adjacency=adj_use, use_weights=True, random_state=SEED
        )
        ari = scib.metrics.ari(self.adata, cluster_key=key_added, label_key="celltype")
        nmi = scib.metrics.nmi(self.adata, cluster_key=key_added, label_key="celltype")
        return ari, nmi

# --- 4. Execution (Baseline vs TopoLa 0.01) ---
topo = TopoLaHandler(adata)
FIXED_LAMBDA = 0.01

# 1. Baseline
ari_base, nmi_base = topo.run_louvain(res=1.0, key_added="louvain_base")

# 2. TopoLa (Fixed Lambda)
ari_topo, nmi_topo = topo.run_louvain(res=1.0, key_added="louvain_topo", alpha=FIXED_LAMBDA)

# --- 5. Results & Plot ---
print(f"\n=== Results (Lambda = {FIXED_LAMBDA}) ===")
print(f"Baseline:       ARI={ari_base:.4f}, NMI={nmi_base:.4f}")
print(f"TopoLa:  ARI={ari_topo:.4f}, NMI={nmi_topo:.4f}")
print(f"Delta:          ARI={ari_topo-ari_base:+.4f}, NMI={nmi_topo-nmi_base:+.4f}")

sc.pl.umap(
    adata, 
    color=["str_batch", "celltype", "louvain_base", "louvain_topo"], 
    legend_loc="on data", 
    legend_fontoutline=2,
    ncols=2,
    title=[
        "Batch", "GT Celltype", 
        f"Baseline (ARI={ari_base:.3f})", 
        f"TopoLa λ={FIXED_LAMBDA} (ARI={ari_topo:.3f})"
    ]
)