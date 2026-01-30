from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any, Literal
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from natsort import natsorted
from packaging import version
from scipy.sparse import csr_matrix
from scanpy import _utils
from scanpy import logging as logg
from scanpy._compat import old_positionals
from scanpy._utils import _choose_graph
from scanpy.tools._utils_clustering import rename_groups, restrict_adjacency

if TYPE_CHECKING:
    from anndata import AnnData
    from scipy.sparse import spmatrix

try:
    from louvain.VertexPartition import MutableVertexPartition
except ImportError:

    class MutableVertexPartition:
        pass

    MutableVertexPartition.__module__ = "louvain.VertexPartition"


@old_positionals(
    "random_state",
    "restrict_to",
    "key_added",
    "adjacency",
    "flavor",
    "directed",
    "use_weights",
    "partition_type",
    "partition_kwargs",
    "neighbors_key",
    "obsp",
    "copy",
)
def louvain(
    adata: AnnData,
    resolution: float | None = None,
    *,
    random_state: _utils.AnyRandom = 0,
    restrict_to: tuple[str, Sequence[str]] | None = None,
    key_added: str = "louvain",
    adjacency: spmatrix | None = None,
    flavor: Literal["vtraag", "igraph", "rapids"] = "vtraag",
    directed: bool = True,
    use_weights: bool = False,
    partition_type: type[MutableVertexPartition] | None = None,
    partition_kwargs: Mapping[str, Any] = {},
    neighbors_key: str | None = None,
    obsp: str | None = None,
    copy: bool = False,
    TopoLa_bool: bool = False,
    TopoLa_alpha: float | None = None,
) -> AnnData | None:
    """... (original docstring unchanged) ..."""
    partition_kwargs = dict(partition_kwargs)
    start = logg.info("running Louvain clustering")
    if (flavor != "vtraag") and (partition_type is not None):
        raise ValueError(
            "partition_type is only a valid argument " 'when flavour is "vtraag"'
        )
    adata = adata.copy() if copy else adata
    if adjacency is None:
        adjacency = _choose_graph(adata, obsp, neighbors_key)

    if restrict_to is not None:
        restrict_key, restrict_categories = restrict_to
        adjacency, restrict_indices = restrict_adjacency(
            adata,
            restrict_key,
            restrict_categories=restrict_categories,
            adjacency=adjacency,
        )
    if flavor in {"vtraag", "igraph"}:
        if flavor == "igraph" and resolution is not None:
            logg.warning('resolution parameter has no effect for flavor "igraph"')
        if directed and flavor == "igraph":
            directed = False
        if not directed:
            logg.debug("    using the undirected graph")

        adjacency_dense = adjacency.toarray()
        adjacency_sparse_df = pd.DataFrame(adjacency_dense)
        adjacency_sparse_df.to_csv("adjacency_sparse.csv", index=False)
        adata.obsm["X_scGPT"] = adjacency_dense
        if TopoLa_bool:  
            adjacency, weights,Matrix_normalized = TopoLa(adjacency, TopoLa_alpha)
            adata.obsm["X_scGPT"] = Matrix_normalized
        else:
            weights = None

        g = _utils.get_igraph_from_adjacency(adjacency, directed=directed)

        if  weights is not None:
            partition_kwargs["weights"] = weights
        else:
            weights = None

        if flavor == "vtraag":
            import louvain

            if partition_type is None:
                partition_type = louvain.RBConfigurationVertexPartition
            if resolution is not None:
                partition_kwargs["resolution_parameter"] = resolution
            if use_weights:
                partition_kwargs["weights"] = weights
            if version.parse(louvain.__version__) < version.parse("0.7.0"):
                louvain.set_rng_seed(random_state)
            else:
                partition_kwargs["seed"] = random_state
            logg.info('    using the "louvain" package of Traag (2017)')
            part = louvain.find_partition(
                g,
                partition_type,
                **partition_kwargs,
            )
            # adata.uns['louvain_quality'] = part.quality()
        else:
            part = g.community_multilevel(weights=weights)
        groups = np.array(part.membership)
    elif flavor == "rapids":
        msg = (
            "flavor='rapids' is deprecated. "
            "Use rapids_singlecell.tl.louvain instead."
        )
        warnings.warn(msg, FutureWarning)
        # nvLouvain only works with undirected graphs,
        # and adjacency must have a directed edge in both directions
        import cudf
        import cugraph

        offsets = cudf.Series(adjacency.indptr)
        indices = cudf.Series(adjacency.indices)
        if use_weights:
            sources, targets = adjacency.nonzero()
            weights = adjacency[sources, targets]
            if isinstance(weights, np.matrix):
                weights = weights.A1
            weights = cudf.Series(weights)
        else:
            weights = None
        g = cugraph.Graph()

        if hasattr(g, "add_adj_list"):
            g.add_adj_list(offsets, indices, weights)
        else:
            g.from_cudf_adjlist(offsets, indices, weights)

        logg.info('    using the "louvain" package of rapids')
        if resolution is not None:
            louvain_parts, _ = cugraph.louvain(g, resolution=resolution)
        else:
            louvain_parts, _ = cugraph.louvain(g)
        groups = (
            louvain_parts.to_pandas()
            .sort_values("vertex")[["partition"]]
            .to_numpy()
            .ravel()
        )
    elif flavor == "taynaud":
        # this is deprecated
        import community

        g = nx.Graph(adjacency)
        partition = community.best_partition(g)
        groups = np.zeros(len(partition), dtype=int)
        for k, v in partition.items():
            groups[k] = v
    else:
        raise ValueError('flavor needs to be "vtraag" or "igraph" or "taynaud".')
    if restrict_to is not None:
        if key_added == "louvain":
            key_added += "_R"
        groups = rename_groups(
            adata,
            key_added=key_added,
            restrict_key=restrict_key,
            restrict_categories=restrict_categories,
            restrict_indices=restrict_indices,
            groups=groups,
        )
    adata.obs[key_added] = pd.Categorical(
        values=groups.astype("U"),
        categories=natsorted(map(str, np.unique(groups))),
    )
    adata.uns[key_added] = {}
    adata.uns[key_added]["params"] = dict(
        resolution=resolution,
        random_state=random_state,
    )
    logg.info(
        "    finished",
        time=start,
        deep=(
            f"found {len(np.unique(groups))} clusters and added\n"
            f"    {key_added!r}, the cluster labels (adata.obs, categorical)"
        ),
    )
    return adata if copy else None

# Sample usage:
# adata = ...  # Load or create your AnnData object
# louvain(adata, TopoLa_bool=True, TopoLa_alpha=0.1)

def TopoLa(A, lambda_val):
    A = A.toarray()
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    S_ = np.diag(S)

    for i in range(len(S)):
        if S[i] != 0:
            S_[i, i] = S[i] ** 3 / (S[i] ** 2 + 1/lambda_val)

    Matrix = np.dot(U, np.dot(S_, Vt))
    
    Matrix_min = Matrix.min()
    Matrix_max = Matrix.max()
    Matrix_normalized = (Matrix - Matrix_min) / (Matrix_max - Matrix_min)
    Matrix_normalized[Matrix_normalized < 1e-3] = 0

    # Save Matrix_normalized to CSV in the current directory
    matrix_normalized_df = pd.DataFrame(Matrix_normalized)
    matrix_normalized_df.to_csv("TopoLa.csv", index=False)
    
    # Generate binary adjacency matrix
    adjacency = (Matrix_normalized != 0).astype(int)
    
    # Extract the weights of non-zero elements
    weights = Matrix_normalized[Matrix_normalized != 0]

    # Convert the adjacency matrix back to sparse format
    adjacency_sparse = csr_matrix(adjacency)
    
    # Convert sparse matrix to dense for CSV export
    adjacency_dense = adjacency_sparse.toarray()
    
    # Save adjacency_sparse (as dense matrix) to CSV in the current directory
    
    
    return adjacency, weights, Matrix_normalized

