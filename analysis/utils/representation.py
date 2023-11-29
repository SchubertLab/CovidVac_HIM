import scanpy as sc
import pandas as pd


tcr_gene_prefixs = ['TRAV', 'TRAJ', 'TRAC', 'TRB', 'TRDV', 'TRDC', 'TRG']


def filter_tcr_genes(adata):
    non_tcr_genes = adata.var_names
    for prefix in tcr_gene_prefixs:
        non_tcr_genes = [el for el in non_tcr_genes if not el.startswith(prefix)]
    adata = adata[:, non_tcr_genes]
    return adata


def filter_high_var(adata, n):
    sc.pp.highly_variable_genes(adata, n_top_genes=n)
    adata = adata[:, adata.var['highly_variable']].copy()
    return adata


def calculate_umap(adata, n_high_var=5000, remove_tcr_genes=True):
    """
    Calculates UMAP coordinates
    :param adata: AnnData object
    :param n_high_var: how many variable genes to use, None uses raw data
    :param remove_tcr_genes: whether to remove TCR genes
    :return: results stored in adata object
    """
    adata_tmp = adata.copy()

    if remove_tcr_genes:
        adata_tmp = filter_tcr_genes(adata_tmp)
    if n_high_var is not None:
        adata_tmp = filter_high_var(adata_tmp, n_high_var)

    sc.pp.neighbors(adata_tmp)
    sc.tl.umap(adata_tmp)
    adata.obsm['X_umap'] = adata_tmp.obsm['X_umap']


def calculate_leiden(adata, resolution, n_high_var=5000, remove_tcr_genes=True):
    adata_tmp = adata.copy()

    if remove_tcr_genes:
        adata_tmp = filter_tcr_genes(adata_tmp)
    if n_high_var is not None:
        adata_tmp = filter_high_var(adata_tmp, n_high_var)

    sc.pp.neighbors(adata_tmp)
    sc.tl.leiden(adata_tmp, resolution=resolution)
    adata.obs['leiden'] = adata_tmp.obs['leiden']


def calculate_umap_by_condition(adata, col_condition, n_high_var=5000, remove_tcr_genes=True):
    df_coords = []
    for s in adata.obs[col_condition].unique():
        adata_tmp = adata[adata.obs[col_condition] == s].copy()
        calculate_umap(adata_tmp, n_high_var, remove_tcr_genes)

        df_tmp = pd.DataFrame(data=adata_tmp.obsm['X_umap'], index=adata_tmp.obs.index, columns=['UMAP1', 'UMAP2'])
        df_coords.append(df_tmp)
    df_coords = pd.concat(df_coords)
    df_coords = df_coords.loc[adata.obs.index]
    adata.obsm[f'X_umap_per_{col_condition}'] = df_coords.values


def calculate_diffmap(adata, n_high_var=5000, remove_tcr_genes=True):
    adata_tmp = adata.copy()

    if remove_tcr_genes:
        adata_tmp = filter_tcr_genes(adata_tmp)
    if n_high_var is not None:
        adata_tmp = filter_high_var(adata_tmp, n_high_var)

    sc.pp.neighbors(adata_tmp)
    sc.tl.diffmap(adata_tmp)
    adata.obsm['X_diffmap'] = adata_tmp.obsm['X_diffmap']
    adata.uns['diffmap_evals'] = adata_tmp.uns['diffmap_evals']


def calculate_dpt(adata, n_branchings=0, n_high_var=5000, remove_tcr_genes=True):
    adata_tmp = adata.copy()

    if remove_tcr_genes:
        adata_tmp = filter_tcr_genes(adata_tmp)
    if n_high_var is not None:
        adata_tmp = filter_high_var(adata_tmp, n_high_var)

    sc.pp.neighbors(adata_tmp)
    sc.tl.dpt(adata_tmp, n_branchings=n_branchings)
    adata.obs['dpt_pseudotime'] = adata_tmp.obs['dpt_pseudotime']
