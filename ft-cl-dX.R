# %%

R

# %%

library(celda)
library(SingleCellExperiment)
library(Seurat)
library(anndataR)

# %%

setwd("/home/liuzhiyu/Projects/neo_caste/Hsal-data-check")

# %%

Hsal_adata <- read_h5ad("./Hsal_ft-cl.h5ad")

# %%

Anndata_2_SCE <- function(anndata) {
    counts_mtx <- Matrix::t(anndata$layers[["counts"]]) # 根据pbmc的例子, decontX需要输入原始计数矩阵(未经过normalize)
    
    if (!inherits(counts_mtx, "dgCMatrix")) {
        counts_mtx <- as(counts_mtx, "CsparseMatrix")
    }
    
    sce <- SingleCellExperiment(
        assays = list(counts = counts_mtx),
        colData = as.data.frame(anndata$obs),
        rowData = as.data.frame(anndata$var)
    )

    rownames(sce) <- rownames(anndata$var)
    colnames(sce) <- rownames(anndata$obs)
    
    return(sce)
}

# %%

Hsal_sce <- Anndata_2_SCE(Hsal_adata)

# %%

Hsal_sce_dX <- decontX(x = Hsal_sce, z = Hsal_sce$leiden_0.5)

# %%

summary(Hsal_sce_dX$decontX_contamination)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0005449 0.2186918 0.3239700 0.3763215 0.5122580 0.9877508 

# %%

materialize <- function(sce, assay_name) {
    # Matrix::t(as(assay(sce, assay_name), "CsparseMatrix"))
    mtx <- assay(sce, assay_name)
    if (!inherits(mtx, "dgCMatrix")) {
        mtx <- as(mtx, "CsparseMatrix")
    }
    return(Matrix::t(mtx))
}

SCE_2_Anndata <- function(sce) {
    counts_t <- materialize(sce, "counts")
    decontXcounts_t <- materialize(sce, "decontXcounts")
    decontX_meta <- metadata(sce)$decontX
    decontX_meta$runParams$logfile <- NULL
    # decontX_meta$runParams$z <- NULL
    anndata <- AnnData(
        X      = counts_t,
        obs    = as.data.frame(colData(sce)),
        var    = as.data.frame(rowData(sce)),
        layers = list(decontXcounts = decontXcounts_t),
        obsm   = list(X_decontX_UMAP = reducedDim(sce, "decontX_UMAP")),
        uns    = list(decontX = decontX_meta)
    )
    return(anndata)
}

# %%

Hsal_sce_dX_adata <- SCE_2_Anndata(Hsal_sce_dX)

# %%

write_h5ad(Hsal_sce_dX_adata, "Hsal_ft-cl-dX_F1_0.5.h5ad", mode = "w")
