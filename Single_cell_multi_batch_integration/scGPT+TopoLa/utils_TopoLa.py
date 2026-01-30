import functools
import json
import logging
import os
from pathlib import Path
import random
import subprocess
from typing import Dict, List, Mapping, Optional, Tuple, Union

import numpy as np
import torch
import pandas as pd
from anndata import AnnData
import scib
from matplotlib import pyplot as plt
from matplotlib import axes
from IPython import get_ipython

import scipy.io
from scib.utils import check_adata, check_batch
from scib.metrics.ari import ari
from scib.metrics.cell_cycle import cell_cycle
from clustring_TopoLa import opt_louvain
from scib.metrics.graph_connectivity import graph_connectivity
from scib.metrics.highly_variable_genes import hvg_overlap
from scib.metrics.isolated_labels import isolated_labels
from scib.metrics.kbet import kBET
from scib.metrics.lisi import clisi_graph, ilisi_graph
from scib.metrics.nmi import nmi
from scib.metrics.pcr import pcr_comparison
from scib.metrics.silhouette import silhouette, silhouette_batch
from scib.metrics.trajectory import trajectory_conservation

from scgpt import logger

def eval_scib_metrics(
    adata: AnnData,
    batch_key: str = "str_batch",
    label_key: str = "celltype",
    notes: Optional[str] = None,
    TopoLa_bool: bool = False,
    TopoLa_alpha: float = 1e-2,
) -> Dict:
    results = metrics(
        adata,
        adata_int=adata,
        batch_key=batch_key,
        label_key=label_key,
        embed="X_scGPT",
        isolated_labels_asw_=False,
        silhouette_=False,
        hvg_score_=False,
        graph_conn_=False,
        pcr_=False,
        isolated_labels_f1_=False,
        trajectory_=False,
        nmi_=True,  # 使用聚类，偏向最佳匹配
        ari_=True,  # 使用聚类，偏向最佳匹配
        cell_cycle_=False,
        kBET_=False,  # kBET 有时返回 nan，需要检查
        ilisi_=False,
        clisi_=False,
        TopoLa_bool=TopoLa_bool,
        TopoLa_alpha=TopoLa_alpha,
    )

    if notes is not None:
        logger.info(f"{notes}")

    logger.info(f"{results}")

    result_dict = results[0].to_dict()
    logger.info(
        "Biological Conservation Metrics: \n"
        f"ASW (cell-type): {result_dict['ASW_label']:.4f}, graph cLISI: {result_dict['cLISI']:.4f}, "
        f"isolated label silhouette: {result_dict['isolated_label_silhouette']:.4f}, \n"
        "Batch Effect Removal Metrics: \n"
        f"PCR_batch: {result_dict['PCR_batch']:.4f}, ASW (batch): {result_dict['ASW_label/batch']:.4f}, "
        f"graph connectivity: {result_dict['graph_conn']:.4f}, graph iLISI: {result_dict['iLISI']:.4f}"
    )

    result_dict["avg_bio"] = np.mean(
        [
            result_dict["NMI_cluster/label"],
            result_dict["ARI_cluster/label"],
            result_dict["ASW_label"],
        ]
    )

    # 移除 result_dict 中的 nan 值
    result_dict = {k: v for k, v in result_dict.items() if not np.isnan(v)}

    return result_dict


def metrics(
    adata,
    adata_int,
    batch_key,
    label_key,
    embed="X_scGPT",
    cluster_key="cluster",
    cluster_nmi=None,
    ari_=False,
    nmi_=False,
    nmi_method="arithmetic",
    nmi_dir=None,
    silhouette_=False,
    si_metric="euclidean",
    pcr_=False,
    cell_cycle_=False,
    organism="mouse",
    hvg_score_=False,
    isolated_labels_=False,  # backwards compatibility
    isolated_labels_f1_=False,
    isolated_labels_asw_=False,
    n_isolated=None,
    graph_conn_=False,
    trajectory_=False,
    kBET_=False,
    lisi_graph_=False,
    ilisi_=False,
    clisi_=False,
    subsample=0.5,
    n_cores=1,
    type_=None,
    verbose=False,
    TopoLa_bool=False,
    TopoLa_alpha=1e-2,
):
    """Master metrics function

    Wrapper for all metrics used in the study.
    Compute of all metrics given unintegrated and integrated anndata object

    :param adata:
        unintegrated, preprocessed anndata object
    :param adata_int:
        integrated anndata object
    :param batch_key:
        name of batch column in adata.obs and adata_int.obs
    :param label_key:
        name of biological label (cell type) column in adata.obs and adata_int.obs
    :param embed:
        embedding representation of adata_int

        Used for:

            + silhouette scores (label ASW, batch ASW),
            + PC regression,
            + cell cycle conservation,
            + isolated label scores, and
            + kBET
    :param cluster_key:
        name of column to store cluster assignments. Will be overwritten if it exists
    :param cluster_nmi:
        Where to save cluster resolutions and NMI for optimal clustering
        If None, these results will not be saved
    :param `ari_`:
        whether to compute ARI using :func:`~scib.metrics.ari`
    :param `nmi_`:
        whether to compute NMI using :func:`~scib.metrics.nmi`
    :param nmi_method:
        which implementation of NMI to use
    :param nmi_dir:
        directory of NMI code for some implementations of NMI
    :param `silhouette_`:
        whether to compute the average silhouette width scores for labels and batch
        using :func:`~scib.metrics.silhouette` and :func:`~scib.metrics.silhouette_batch`
    :param si_metric:
        which distance metric to use for silhouette scores
    :param `pcr_`:
        whether to compute principal component regression using :func:`~scib.metrics.pc_comparison`
    :param `cell_cycle_`:
        whether to compute cell cycle score conservation using :func:`~scib.metrics.cell_cycle`
    :param organism:
        organism of the datasets, used for computing cell cycle scores on gene names
    :param `hvg_score_`:
        whether to compute highly variable gene conservation using :func:`~scib.metrics.hvg_overlap`
    :param `isolated_labels_`:
        whether to compute both isolated label scores using :func:`~scib.metrics.isolated_labels`
    :param `isolated_labels_f1_`:
        whether to compute isolated label score based on F1 score of clusters vs labels using
        :func:`~scib.metrics.isolated_labels`
    :param `isolated_labels_asw_`:
        whether to compute isolated label score based on ASW (average silhouette width) using
        :func:`~scib.metrics.isolated_labels`
    :param `n_isolated`:
        maximum number of batches per label for label to be considered as isolated
    :param `graph_conn_`:
        whether to compute graph connectivity score using :func:`~scib.metrics.graph_connectivity`
    :param `trajectory_`:
        whether to compute trajectory score using :func:`~scib.metrics.trajectory_conservation`
    :param `kBET_`:
        whether to compute kBET score using :func:`~scib.metrics.kBET`
    :param `lisi_graph_`:
        whether to compute both cLISI and iLISI using :func:`~scib.metrics.lisi_graph`
    :param `clisi_`:
        whether to compute cLISI using :func:`~scib.metrics.clisi_graph`
    :param `ilisi_`:
        whether to compute iLISI using :func:`~scib.metrics.ilisi_graph`
    :param subsample:
        subsample fraction for LISI scores
    :param n_cores: number of cores to be used for LISI functions
    :param `type_`:
        one of 'full', 'embed' or 'knn' (used for kBET and LISI scores)
    """

    check_adata(adata)
    check_batch(batch_key, adata.obs)
    check_batch(label_key, adata.obs)

    check_adata(adata_int)
    check_batch(batch_key, adata_int.obs)
    check_batch(label_key, adata_int.obs)
    scipy.io.savemat("labs.mat", {label_key:label_key})
    # clustering
    
    if nmi_ or ari_:
        res_max, nmi_max, nmi_all = opt_louvain(
            adata_int,
            label_key=label_key,
            cluster_key=cluster_key,
            use_rep=embed,
            function=nmi,
            plot=False,
            verbose=verbose,
            inplace=True,
            force=True,
            TopoLa_bool=TopoLa_bool,
            TopoLa_alpha=TopoLa_alpha,
        )
        if cluster_nmi is not None:
            nmi_all.to_csv(cluster_nmi, header=False)
            print(f"saved clustering NMI values to {cluster_nmi}")

    if nmi_:
        print("NMI...")
        nmi_score = nmi(
            adata_int,
            cluster_key=cluster_key,
            label_key=label_key,
            implementation=nmi_method,
            nmi_dir=nmi_dir,
        )
    else:
        nmi_score = np.nan

    if ari_:
        print("ARI...")
        ari_score = ari(adata_int, cluster_key=cluster_key, label_key=label_key)
    else:
        ari_score = np.nan

    if silhouette_:
        print("Silhouette score...")
        # global silhouette coefficient
        asw_label = silhouette(
            adata_int, label_key=label_key, embed=embed, metric=si_metric
        )
        # silhouette coefficient per batch
        asw_batch = silhouette_batch(
            adata_int,
            batch_key=batch_key,
            label_key=label_key,
            embed=embed,
            metric=si_metric,
            return_all=False,
            verbose=False,
        )
    else:
        asw_label = np.nan
        asw_batch = np.nan

    if pcr_:
        print("PC regression...")
        pcr_score = pcr_comparison(
            adata, adata_int, embed=embed, covariate=batch_key, verbose=verbose
        )
    else:
        pcr_score = np.nan

    if cell_cycle_:
        print("cell cycle effect...")
        cc_score = cell_cycle(
            adata,
            adata_int,
            batch_key=batch_key,
            embed=embed,
            agg_func=np.mean,
            organism=organism,
        )
    else:
        cc_score = np.nan

    if isolated_labels_f1_ or isolated_labels_:
        print("Isolated labels F1...")
        il_score_f1 = isolated_labels(
            adata_int,
            label_key=label_key,
            batch_key=batch_key,
            embed=embed,
            cluster=True,
            iso_threshold=n_isolated,
            verbose=False,
        )
    else:
        il_score_f1 = np.nan

    if isolated_labels_asw_ or isolated_labels_:
        print("Isolated labels ASW...")
        il_score_asw = (
            isolated_labels(
                adata_int,
                label_key=label_key,
                batch_key=batch_key,
                embed=embed,
                cluster=False,
                iso_threshold=n_isolated,
                verbose=False,
            )
            if silhouette_
            else np.nan
        )
    else:
        il_score_asw = np.nan

    if graph_conn_:
        print("Graph connectivity...")
        graph_conn_score = graph_connectivity(adata_int, label_key=label_key)
    else:
        graph_conn_score = np.nan

    if kBET_:
        print("kBET...")
        kbet_score = kBET(
            adata_int,
            batch_key=batch_key,
            label_key=label_key,
            type_=type_,
            embed=embed,
            scaled=True,
            verbose=verbose,
        )
    else:
        kbet_score = np.nan

    if lisi_graph_:
        clisi_ = True
        ilisi_ = True

    if clisi_:
        print("cLISI score...")
        clisi = clisi_graph(
            adata_int,
            batch_key=batch_key,
            label_key=label_key,
            type_=type_,
            subsample=subsample * 100,
            scale=True,
            n_cores=n_cores,
            verbose=verbose,
        )
    else:
        clisi = np.nan

    if ilisi_:
        print("iLISI score...")
        ilisi = ilisi_graph(
            adata_int,
            batch_key=batch_key,
            type_=type_,
            subsample=subsample * 100,
            scale=True,
            n_cores=n_cores,
            verbose=verbose,
        )
    else:
        ilisi = np.nan

    if hvg_score_:
        hvg_score = hvg_overlap(adata, adata_int, batch_key)
    else:
        hvg_score = np.nan

    if trajectory_:
        print("Trajectory conservation score...")
        trajectory_score = trajectory_conservation(
            adata,
            adata_int,
            label_key=label_key,
            # batch_key=batch_key
        )
    else:
        trajectory_score = np.nan

    results = {
        "NMI_cluster/label": nmi_score,
        "ARI_cluster/label": ari_score,
        "ASW_label": asw_label,
        "ASW_label/batch": asw_batch,
        "PCR_batch": pcr_score,
        "cell_cycle_conservation": cc_score,
        "isolated_label_F1": il_score_f1,
        "isolated_label_silhouette": il_score_asw,
        "graph_conn": graph_conn_score,
        "kBET": kbet_score,
        "iLISI": ilisi,
        "cLISI": clisi,
        "hvg_overlap": hvg_score,
        "trajectory": trajectory_score,
    }

    return pd.DataFrame.from_dict(results, orient="index")
