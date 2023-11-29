import scanpy as sc

import seaborn as sb
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import os
import sys
sys.path.append(os.path.dirname(__file__))

import annotation as anno


def separate_umaps_by_condition(adata, condition, n_rows, n_cols, do_int_sort=False, title=None, **kwargs):
    classes = adata.obs[condition].unique().tolist()
    for nan in [None, np.nan, 'nan']:
        if nan in classes:
            classes.remove(nan)
    classes = sorted(classes)
    if do_int_sort:
        classes = sorted([int(el) for el in classes])
        classes = [str(el) for el in classes]
    assert len(classes) <= n_rows * n_cols, f'n_rows*n_cols={n_rows*n_cols} does not fit {len(classes)} classes.'

    fig, axes = plt.subplots(ncols=n_cols, nrows=n_rows, figsize=(n_cols * 3, n_rows * 3))
    axes = axes.reshape(-1)

    for i, cluster in enumerate(classes):
        sc.pl.umap(adata, color=condition, groups=cluster, ax=axes[i], show=False,
                   palette=['tab:orange'] * len(classes),
                   title=cluster, **kwargs)
        axes[i].get_legend().remove()
        axes[i].set_ylabel(None)
        axes[i].set_xlabel(None)

    for j in range(i+1, n_rows*n_cols):
        axes[j].axis('off')
    fig.tight_layout()
    plt.suptitle(title)
    plt.show()
    adata.uns.pop(f'{condition}_colors')


def distributions_over_columns(adata, columns, n_rows, n_cols, y_lim=0.003, x_lim=None, title=None):
    assert len(columns) <= n_rows * n_cols, f'n_rows*n_cols={n_rows*n_cols} does not fit {len(columns)} classes.'
    fig, axes = plt.subplots(ncols=n_cols, nrows=n_rows, figsize=(n_cols * 3, n_rows * 3))
    if type(axes) == np.ndarray:
        axes = axes.reshape(-1)
    else:
        axes = [axes]

    for i, ht in enumerate(columns):
        if np.sum(np.isnan(adata.obs[ht])) != len(adata):
            sb.distplot(adata[adata.obs[ht] >= 0].obs[ht], ax=axes[i])
        axes[i].set_ylim([0, y_lim])
        if x_lim is not None and x_lim < max(adata.obs[ht]):
            axes[i].set_xlim([0, x_lim])
        axes[i].set_title(ht)
    fig.suptitle(title)
    fig.tight_layout()
    plt.show()


def plot_marker_genes(adata):
    for name, genes in anno.marker_collection.items():
        genes_present = [el for el in genes if el in adata.var_names] + [name]
        sc.pl.umap(adata, color=genes_present)
        plt.show()


def adt_counts_by_condition(adata, collumns_counts, collumns_label,  n_rows, n_cols, do_log=True, vmax=None):
    fig = plt.figure(figsize=(n_cols * 3, n_rows * 3))

    gs = fig.add_gridspec(figure=fig, nrows=n_rows, ncols=n_cols, hspace=0.5, )

    for i, c in enumerate(collumns_counts):
        gs_internal = gridspec.GridSpecFromSubplotSpec(subplot_spec=gs[i], nrows=2, ncols=1)
        axes_0 = fig.add_subplot(gs_internal[0])
        axes_1 = fig.add_subplot(gs_internal[1])

        vmax = max(adata.obs[c].max() + 20, 10) if vmax is None else vmax

        if np.sum(adata.obs[collumns_label] == c) > 0:
            sb.distplot(adata[adata.obs[collumns_label] == c].obs[c], ax=axes_0, color='tab:pink',
                        hist=True, kde_kws={'fill': True, 'bw_adjust': 1}, axlabel=False, kde=True)
            axes_0.set_title(None)
            if do_log:
                axes_0.set_xscale('symlog')
            axes_0.set_xlim((0, vmax))
            axes_0.set_xticks([])
            axes_0.set_yticks([])
            axes_0.set_ylabel(None)
        else:
            axes_0.set_visible(False)  # axis('off')
            pass

        sb.distplot(adata[adata.obs[collumns_label] != c].obs[c], ax=axes_1, color='tab:olive',
                    hist=True, kde_kws={'fill': True}, kde=True)
        axes_1.set_yticks([])
        axes_1.set_ylabel(None)
        # axes[1].set_yticklabels('')
        if do_log:
            axes_1.set_xscale('symlog')
        axes_1.set_xlim((0, vmax))

    plt.tight_layout()
    plt.show()


def plot_largest_clonotypes(adata, title, n, ncols, nrows):
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(ncols*3, nrows*3))
    axes = axes.reshape(-1)

    adata.obs['clone_id'] = adata.obs['clone_id'].astype(str)
    top_35_clonotypes = list(adata.obs['clone_id'].value_counts().head(n+1).index)
    if 'nan' in top_35_clonotypes:
        top_35_clonotypes.remove('nan')
    top_35_clonotypes = top_35_clonotypes[:n]

    adata_tmp = adata.copy()
    for cluster in top_35_clonotypes:
        adata_tmp.obs[f'clonotype_{cluster}'] = (adata_tmp.obs['clone_id'] == cluster).astype(str)

    for i, cluster in enumerate(top_35_clonotypes):
        sc.pl.umap(adata_tmp, color=f'clonotype_{cluster}', groups='True', ax=axes[i], show=False,
                   palette=['tab:red'] * n, title=cluster)
        axes[i].get_legend().remove()

    fig.suptitle(title)
    fig.tight_layout()
    plt.show()


def umap_with_background_condition(adata, color, col_condition, condition, axes, size=None,
                                   highlight=False, factor=3, adata_back=None, edgecolor='black', **kwargs):
    if adata_back is None:
        sc.pl.umap(adata, palette=['silver'], size=size, ax=axes, show=False)  # todo: shift cmap + nacolor to utils
    else:
        sc.pl.umap(adata_back, palette=['silver'], size=size, ax=axes, show=False)

    adata_tmp = adata[adata.obs[col_condition] == condition]
    if not highlight:
        sc.pl.umap(adata_tmp, color=color, ax=axes, show=False, **kwargs)
    else:
        sc.pl.umap(adata_tmp[adata_tmp.obs[highlight] != 'True'], color=color, ax=axes, show=False, size=size, **kwargs)
        cb = axes.collections[-1].colorbar
        if cb is not None:
            cb.remove()
        sc.pl.umap(adata_tmp[adata_tmp.obs[highlight] == 'True'], color=color, ax=axes, show=False, size=factor*size,
                   edgecolor=edgecolor, linewidths=0.5,
                   **kwargs)
