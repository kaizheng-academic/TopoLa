import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import shutil
from GraphST import GraphST
from src.utils_TopoLa import clustering

# Run device, by default, the package is implemented on 'cpu'. We recommend using GPU.
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

# List of datasets and their respective number of clusters
datasets = {
    '151507': 7, '151508': 7, '151509': 7, '151510': 7,
    '151669': 5, '151670': 5, '151671': 5, '151672': 5,
    '151673': 7, '151674': 7, '151675': 7, '151676': 7
}

# Fixed lambda value
lambda_val = 0.1

# Ensure embedding directory exists
embedding_dir = './embedding/'
os.makedirs(embedding_dir, exist_ok=True)

def run_clustering_and_save_results():
    """
    Perform clustering on each dataset with a fixed lambda value and save the embeddings, ARI, and NMI results.
    """
    # Create a DataFrame to store ARI and NMI results
    results = pd.DataFrame(columns=['Dataset', 'ARI', 'NMI'])

    for dataset, n_clusters in datasets.items():
        print(f"Processing dataset {dataset} with {n_clusters} clusters")
        # Read data
        file_fold = './data/' + str(dataset)  # Replace with your download path
        adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True)
        adata.var_names_make_unique()
        
        # Define and train model
        model = GraphST.GraphST(adata, device=0)
        adata = model.train()

        # Set radius for neighbors
        radius = 50

        # Add ground_truth
        df_meta = pd.read_csv(file_fold + '/metadata.tsv', sep='\t')
        adata.obs['ground_truth'] = df_meta['layer_guess'].values

        adata = adata[~pd.isnull(adata.obs['ground_truth'])]

        # Clustering with the fixed lambda value
        clustering(adata, n_clusters, radius=radius, method='mclust', refinement=True, lambda_val=lambda_val)

        # Calculate ARI and NMI
        ARI = metrics.adjusted_rand_score(adata.obs['domain'], adata.obs['ground_truth'])
        NMI = metrics.normalized_mutual_info_score(adata.obs['domain'], adata.obs['ground_truth'])
        
        print(f'ARI for dataset {dataset} : {ARI}, NMI: {NMI}')

       

        # Append results to the DataFrame
        results = results.append({'Dataset': dataset, 'ARI': ARI, 'NMI': NMI}, ignore_index=True)

    # Save results to CSV
    results.to_csv("Result.csv", index=False)

    print("Processing complete.")

run_clustering_and_save_results()
