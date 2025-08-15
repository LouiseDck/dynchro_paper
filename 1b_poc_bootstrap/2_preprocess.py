from rpy2.robjects.packages import importr
import scanpy as sc

def assign_iroot(adata):
    adata.uns['iroot'] = np.flatnonzero(adata.obs['milestone'] == 'iroot')[0]
    return adata

def cutoff(adata, milestone):
    adata = adata[adata.obs['milestone'] != milestone]
    return adata

def preprocess(dataset):
    sc.pp.normalize_total(dataset, target_sum=1e4)
    sc.pp.log1p(dataset)
    sc.pp.pca(dataset)
    sc.pp.neighbors(dataset, n_neighbors=10)
    sc.tl.umap(dataset)
    sc.tl.leiden(dataset)
    return dataset

def assign_slingshot(adata, start, ends):
    adata.uns['slingshot_start'] = start
    adata.uns['slingshot_ends'] = ends

    return adata

def select_pseudotime(adata):
    pass

def run_slingshot(adata, start, ends, dimred):
    slingshot = importr("slingshot")
    anndataR = importr("anndataR")
    ss_result = slingshot.slingshot(
        adata.obsm[dimred], 
        cl = adata.obs['leiden'],
        start_clus = adata.uns['slingshot_start'],
        end_clus = adata.uns['slingshot_ends']
    )

    robjects.r["colData"]
    return adata