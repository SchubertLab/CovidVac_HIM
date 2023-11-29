"""
File for creating output tables from the h5ad file.
"""
import pandas as pd


def deg_table(adata, filename, key='rank_genes_groups', prefix=''):
    categories = adata.uns[key]['names'].dtype.names

    dict_deg = {}
    for cluster in categories:
        names = adata.uns[key]['names'][cluster].tolist()
        scores = adata.uns[key]['scores'][cluster].tolist()
        dict_deg[cluster] = list(zip(names, scores))

    df_degs = pd.DataFrame(dict_deg)
    df_degs.columns = [f'{prefix}{el}' for el in df_degs.columns]
    df_degs.to_csv(filename)


def clonotype_by_condition(adata, condition, filepath):
    categories = adata.obs[condition].unique()
    for cat in categories:
        ct_content = adata[(adata.obs[condition] == cat) & ~adata.obs['clone_id'].isna()].obs
        ct_content = ct_content[ct_content['clone_id'] != 'nan']
        ct_content = ct_content[[condition, 'clone_id', 'clonotype_sequence', 'v_genes', 'j_genes',
                                 'clone_size', 'celltype']]
        ct_content = ct_content.drop_duplicates()
        ct_content = ct_content.reset_index(drop=True)
        ct_content.to_csv(f'{filepath}/clonotypes_{condition}_{cat}.csv')
