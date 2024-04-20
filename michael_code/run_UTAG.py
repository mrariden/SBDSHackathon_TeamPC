import anndata
from utag import utag
import numpy as np


adata = anndata.read_h5ad('preprocessed_adata.h5ad')

obs = adata.obs[['fov', 'CenterX_global_px', 'CenterY_global_px']]
X = adata.X.toarray()
X = X.astype(np.float64)
var = adata.var

minimal_adata = anndata.AnnData(
    X = X,
    obs = obs,
    var = var
)

minimal_adata.obsm['spatial'] = np.array(minimal_adata.obs[['CenterY_global_px', 'CenterX_global_px']])

from utag import utag

utag_results = utag(
    minimal_adata,
    slide_key="fov",
    max_dist=15,
    normalization_mode='l1_norm',
    apply_clustering=True,
    clustering_method = 'leiden', 
)


utag_results.write('preprocessed_adata-UTAG.h5ad')
