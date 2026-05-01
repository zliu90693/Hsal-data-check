# %%

import scanpy as sc
from scipy.sparse import csr_matrix

# %%

Hsal_adata = sc.read_text("./GSE135513_d30_DGE.txt").T

# %%


