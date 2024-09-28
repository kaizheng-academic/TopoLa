import pickle
import scanpy as sc
import scib
import pandas as pd
import gc

def evaluate_nr_alpha(filename):
    def load_data(filename):
        with open(filename, 'rb') as f:
            adata, label_key, resolution_key = pickle.load(f)
        print(f'Data loaded from {filename}')
        return adata, label_key, resolution_key

    # 调用load_data函数
    adata, label_key, resolution_key = load_data(filename)

    # 定义不同的NR_alpha值和resolution值
    nr_alpha_values = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
    resolution_values = [0.1 * i for i in range(1, 3)]

    # 初始化最大ARI分数和对应的NR_alpha值、resolution值
    max_ari_score = -1
    best_nr_alpha = None
    best_resolution = None

    # 循环测试不同的NR_alpha值和resolution值
    for resolution in resolution_values:
        for alpha in nr_alpha_values:
            try:
                # 执行louvain聚类
                sc.tl.louvain(adata, resolution=resolution, key_added=resolution_key, NR_bool=True, NR_alpha=alpha)

                # 计算ARI分数
                ari_score = scib.metrics.ari(adata, cluster_key=resolution_key, label_key=label_key)

                # 打印当前的ARI分数、NR_alpha值和resolution值
                print(f'Resolution: {resolution}, NR_alpha: {alpha}, ARI Score: {ari_score}')

                # 更新最大ARI分数和对应的NR_alpha值、resolution值
                if ari_score > max_ari_score:
                    max_ari_score = ari_score
                    best_nr_alpha = alpha
                    best_resolution = resolution
            except Exception as e:
                print(f"Error with Resolution: {resolution}, NR_alpha: {alpha}: {e}")
                continue  # 跳过当前循环

    # 输出最大ARI分数和对应的NR_alpha值、resolution值
    print(f'Maximum ARI Score: {max_ari_score}, Corresponding NR_alpha: {best_nr_alpha}, Corresponding Resolution: {best_resolution}')
    return max_ari_score, best_nr_alpha, best_resolution

# 定义四个文件名
filenames = [
    'covid_saved_data.pkl',
    'pbmc_saved_data.pkl',
    'BMMC_saved_data.pkl'
]

# 循环处理每个文件
results = []
for filename in filenames:
    print(f'\nProcessing {filename}...')
    try:
        max_ari_score, best_nr_alpha, best_resolution = evaluate_nr_alpha(filename)
        results.append({'filename': filename, 'max_ari_score': max_ari_score, 'best_nr_alpha': best_nr_alpha, 'best_resolution': best_resolution})
    except Exception as e:
        print(f"Failed to process {filename}: {e}")
    gc.collect()  # 显式调用垃圾回收，释放内存

# 保存结果到CSV文件
results_df = pd.DataFrame(results)
results_df.to_csv('results.csv', index=False)
