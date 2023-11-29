import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

from tcrdist.repertoire import TCRrep


import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from matplotlib.patches import Patch



def data_preparation(adata, cols=None):
    """
    Transforms the adata into a df of the TCRdist format
    :param adata: input data
    :param cols: additional columns to copy
    :return: df in TCRdist format
    """
    if cols is None:
        cols = []
    df_tcrdist_format = adata.obs[['clone_id', 'IR_VJ_1_v_call', 'IR_VDJ_1_v_call',
                                   'IR_VJ_1_junction_aa', 'IR_VDJ_1_junction_aa'] + cols]
    df_tcrdist_format = df_tcrdist_format.drop_duplicates()

    rename_tcrdist = {
        'IR_VJ_1_v_call': 'v_a_gene',
        'IR_VDJ_1_v_call': 'v_b_gene',
        'IR_VJ_1_junction_aa': 'cdr3_a_aa',
        'IR_VDJ_1_junction_aa': 'cdr3_b_aa',
    }
    df_tcrdist_format = df_tcrdist_format.rename(columns=rename_tcrdist)

    for col in ['v_a_gene', 'v_b_gene']:
        df_tcrdist_format[col] = df_tcrdist_format[col].astype(str)
        df_tcrdist_format[col] = df_tcrdist_format[col] + '*01'

    df_tcrdist_format['count'] = 1
    df_tcrdist_format.index = df_tcrdist_format['clone_id']
    df_tcrdist_format.index.name = None
    return df_tcrdist_format


def get_tcr_distance(df, chains=None):
    """
    Measure TCR distances between clones of the same dataset
    :param: df: df with filtered data obtained from prepare data function
    :return: df_tcr_distances: df with distances orderd by clone_id
    """

    if chains is None or chains == 'both':
        chain_list = ['alpha', 'beta']
    else:
        chain_list = [chains]

    tr = TCRrep(cell_df=df,
                organism='human',
                chains=chain_list,
                compute_distances=False,
                deduplicate=False,
                db_file='alphabeta_gammadelta_db.tsv')
    tr.compute_distances()

    if chains is None or chains == 'both':
        dist_total = tr.pw_alpha + tr.pw_beta
    elif chains == 'alpha':
        dist_total = tr.pw_alpha
    else:
        dist_total = tr.pw_beta

    index = tr.clone_df['clone_id']
    df_tcr_distances = pd.DataFrame(data=dist_total, index=index, columns=index)

    return df_tcr_distances


def plot_tcrdist_heatmap(df_dists, df_anno=None, colors_cols=None, text_col=None, file_name=None, no_show=False,
                         size_factor=1):
    df_dists = df_dists.copy()
    row_colors = None
    if colors_cols is not None:
        row_colors = {name: [c_dict[df_anno.loc[ct, name]] for ct in df_dists.index]
                      for name, c_dict in colors_cols.items()}
        row_colors = pd.DataFrame(row_colors, index=df_dists.index)

    if text_col is not None:
        df_dists.index = [f'{ct} in {df_anno.loc[ct, text_col]}' for ct in df_dists.index]
        row_colors.index = [f'{ct} in {df_anno.loc[ct, text_col]}' for ct in row_colors.index]

    linkage = hc.linkage(sp.distance.squareform(df_dists), method='average')
    plot = sb.clustermap(df_dists, row_linkage=linkage, col_linkage=linkage,
                         figsize=((40+len(colors_cols)*1)*size_factor, 40*size_factor),
                         row_colors=row_colors, colors_ratio=0.01,
                         xticklabels=1, yticklabels=1)

    if colors_cols is not None:
        labels = [l for c_dicts in colors_cols.values() for l in c_dicts.keys() if l != 'nan']
        handles = [Patch(facecolor=h) for c_dicts in colors_cols.values() for l, h in c_dicts.items() if l != 'nan']
        #        handles = [Patch(facecolor=color_bindings[name]) for name in color_bindings]
        plt.legend(handles, labels, bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
                   fontsize='xx-large')

    if file_name is not None:
        plt.savefig(file_name)
    if no_show:
        return plot
    plt.show()


def dists_to_graph(dists, threshold, percentile=100, min_compenents=1, annotation=None, path_save=None):
    """
    Transform pairwise distances to a graph network
    :param dists: df containing the pairwise distances
    :param threshold: threshold to separate connected components
    :param percentile: percentile of betweenness to be kept
    :param min_compenents: minimal amount of components within a cluster
    :param annotation: df annotating the nodes (index=clone_id, cols=attributes)
    :param path_save: output location for a gefx file
    :return: networkx graph
    """
    import networkx as nx
    edges = dists.copy()
    edges.columns = edges.columns.astype(str)
    edges['clone_id_1'] = edges.index
    edges = edges.melt(id_vars=['clone_id_1'])
    edges = edges.rename(columns={'clone_id_1': 'source',
                                  'clone_id': 'target',
                                  'value': 'weight'})

    edges['target'] = edges['target'].astype(str)
    edges['source'] = edges['source'].astype(str)
    edges = edges[edges['source'] != edges['target']]
    edges = edges[edges['weight'] <= threshold]

    edges['weight'] = 1
    graph = nx.from_pandas_edgelist(edges, edge_attr=True)

    if percentile != 100:
        betweenness = nx.betweenness_centrality(graph)
        percentile_val = np.percentile(list(betweenness.values()), percentile)
        between_nodes = [node for node, val in betweenness.items() if val > percentile_val]
        graph.remove_nodes_from(between_nodes)

    if min_compenents != 1:
        small_clusters = [component for component in list(nx.connected_components(graph))
                          if len(component) < min_compenents]
        small_clusters = [node for component in small_clusters for node in component]
        graph.remove_nodes_from(small_clusters)

    if annotation is not None:
        for col in annotation.columns:
            anno_dict = dict(zip(annotation.index, annotation[col]))
            nx.set_node_attributes(graph, anno_dict, col)

    if path_save:
        nx.write_gexf(graph, path_save)
    return graph
