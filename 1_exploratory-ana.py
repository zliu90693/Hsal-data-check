# %%

import scanpy as sc
from scipy.sparse import csr_matrix

# %%

Hsal_adata = sc.read_text("./GSE135513_d30_DGE.txt").T

# %%

Hsal_adata.var_names_make_unique()

# %%

sc.pp.calculate_qc_metrics(
        Hsal_adata, 
        inplace=True, # 直接修改原对象, 不返回新对象
        percent_top=[20], # 计算表达量排名前20基因占总表达量百分比, 判断表达量是否集中于少数基因
        log1p=True # 对部分指标进行log1p变换
    )

# %%

print(Hsal_adata.obs['total_counts'].describe())

# %%

import numpy as np
import seaborn as sns
from scipy.stats import median_abs_deviation

# %%

def check_3_QC_covariates(adata):
    sns.displot(adata.obs["total_counts"], bins=100, kde=False) 
    sc.pl.violin(adata, 'total_counts')
    # sc.pl.violin(adata, "pct_counts_mt")
    # sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# %%

check_3_QC_covariates(Hsal_adata)

# %%

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

# %%

def add_outlier_column(adata, nmad):
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", nmad) # 注意, 这里的标准是log!
        | is_outlier(adata, "log1p_n_genes_by_counts", nmad)
        | is_outlier(adata, "pct_counts_in_top_20_genes", nmad)
    )
    print(adata.obs["outlier"].value_counts())

# %%

add_outlier_column(Hsal_adata, nmad=5) # Outliner: 12747; Not Outliner: 102641

# %%

def filt_with_outliner_mtoutliner(adata):
    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier)].copy() # 注意, 这里不是原地修改, 这是创建了一个新对象, 如果不返回, 这个对象在函数内会被丢弃
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
    return adata

# %%

Hsal_adata = filt_with_outliner_mtoutliner(Hsal_adata)

# %%

check_3_QC_covariates(Hsal_adata)

# %%

Hsal_adata.write("./Hsal_initalQC.h5ad")

# %%

def neo_pre_clust(adata, reso):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(adata, n_pcs=20)
    
    sc.tl.leiden(adata, key_added=f"leiden_{reso}", resolution=reso)

neo_pre_clust(Hsal_adata, 0.5)

# %%

Hsal_adata.write("./Hsal_ft-cl.h5ad")

# %%

# ft-cl-dX.R 部分

# %%

Hsal_ft_cl_dX = sc.read_h5ad("./Hsal_ft-cl-dX_F1_0.5.h5ad")

# %%

import matplotlib.pyplot as plt

def visualize(adata):
    sc.pl.embedding(adata, basis="decontX_UMAP", color="decontX_contamination", cmap="Reds")
    sc.pl.embedding(adata, basis="decontX_UMAP", color="decontX_clusters")

    plt.hist(adata.obs["decontX_contamination"], bins=50)
    plt.xlabel("Contamination Score")
    plt.ylabel("Cell Count")
    plt.show()

    sc.pl.violin(adata, keys="decontX_contamination", groupby="decontX_clusters")

    print(adata.obs["decontX_contamination"].describe())

# %%

visualize(Hsal_ft_cl_dX)

# %%
