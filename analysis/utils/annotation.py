import numpy as np
import pandas as pd
import scanpy as sc
from collections import Counter
from tqdm import tqdm

tqdm.pandas()


# marker genes
monozytes_classical = ['PGD', 'MS4A6A', 'FCGR2A', 'CSF3R', 'IL17RA', 'DPYD', 'ANPEP', 'CD36', 'ATP6V0E1', 'LMAN2']
monozytes_nonclassical = ['TCF7L2', 'CDKN1C', 'SIGLEC10', 'IER2', 'TOB1', 'NFKBIZ',
                          'NAP1L1', 'FCGR3B', 'CSTA']

nk_markers = ['NCAM1', 'KLRB1', 'IL32', 'S100A4']
markers_nk_b1 = ['FGFBP2', 'GZMB', 'GZMA', 'SPON2', 'S100A4', 'CST7', 'FCGR3A', 'IGFBP7', 'GZMH', 'CFL1']
markers_nk_b2 = ['GZMK', 'CD44', 'PPP1R14B', 'CXCR3', 'RPL36A', 'SCML1', 'COTL1', 'NCF1', 'XCL1', 'HLA-DRB1']

markers_receptor = ['CD3D', 'CD3E', 'TRAC']
markers_gamma_delta = ['TRDC', 'TRGC2', 'TRGC1']

marker_collection = {
    'markers_MonozyteClassical': monozytes_classical,
    'markers_MonozytesNon-Classical': monozytes_nonclassical,
    'markers_NK': nk_markers,
    'markers_NK_b1': markers_nk_b1,
    'markers_NK_b2': markers_nk_b2,
    'markers_TCR': markers_receptor,
    'markers_GammaDeltaTC': markers_gamma_delta,
}


def add_seumois_score(adata):
    """ Add Score described in Seumois et al., Science Immunol. 2020."""
    gene_list = pd.read_csv('../../data/scores/ifn_seumois2020.txt')
    sc.tl.score_genes(adata, gene_list, ctrl_size=100, n_bins=24, score_name='ifn_seumois')


def add_meckiff_score(adata):
    """ Add Score described in Meckiff et al. Cell 2020."""
    score_genes = pd.read_excel('../../data/scores/1-s2.0-S0092867420313076-mmc2.xlsx',
                                sheet_name='H- Gene lists', header=3)

    # many genes in the table are listed with different gene names
    # dictionary to map these genes according to geneCards
    aliases = {
        'FAM46C': 'TENT5C',
        'KIAA0101': 'PCLAF',
        'FAM179A': 'TOGARAM2',
        'TMEM2': 'CEMIP2',
        'GRAMD3': 'GRAMD2B',
        'SEPN1': 'SELENON',
        'FAM134B': 'RETREG1',
        'GPR56': 'ADGRG1',
        'GPR114': 'ADGRG5',
        'BZRAP1': 'TSPOAP1',
        'FAM46A': 'TENT5A',
        'CTB-91J4.1': 'TBC1D3B',
        'PLA2G16': 'PLAAT3',
        'MIR4435-1HG': 'MIR4435-2HG',
        'CTB-58E17.1': 'LOC105371763',
        'IL8': 'CXCL8',
        'MTERFD3': 'MTERF2',
        'C9orf142': 'PAXX',
        'KIAA1715': 'LNPK',
        'LINC00152': 'CYTOR',
        'GATS': 'CASTOR3P',
        'GLTSCR1L': 'BICRAL',
        'C20orf24': 'RAB5IF',
        'C10orf128': 'TMEM273',
        'FAM129A': 'NIBAN1',
        'C4orf32': 'FAM241A',
        'HN1': 'JPT1',
        'MLLT4': 'AFDN',
        'LINC00341': 'SYNE3',
        'PIDD': 'PIDD1',
        'PTPLA': 'HACD1',
        'C15orf53': 'LINC02694',
        'CASC5': 'KNL1',
        'SGOL2': 'SGO1',
        'SGOL1': 'SGO1',
        'ATPIF1': 'ATP5IF1',
        'FLJ22447': 'LINC03033',
        'LOC100131347': 'RDM1P5',
        'WHSC1': 'NSD2'
    }
    score_genes = score_genes.replace(aliases)

    names = ['TFH', 'cytoxicity', 'TH17', 'Type I+II IFN signaling', 'Treg', 'cell cycle', 'Tfrs']
    names = [f'meckiff_{n}' for n in names]
    score_genes.columns = names
    for c in score_genes.columns:
        genes = score_genes[c]
        genes = genes[~genes.isna()].tolist()
        sc.tl.score_genes(adata, genes, ctrl_size=100, n_bins=24, score_name=c)


def add_score(adata, genes, name):
    """
    Adds the average over the normalized gene counts.
    :param adata: AnnData object, contains elements of genes
    :param genes: list of genes
    :param name: column name for storing the scores
    :return: result is stored in adata.obs[name]
    """
    adata.obs[name] = np.mean((
                adata[:, [x for x in genes if x in adata.var_names]].X -
                adata[:, [x for x in genes if x in adata.var_names]].X.mean(axis=0)
            ) / adata[:, [x for x in genes if x in adata.var_names]].X.todense().std(axis=0),
            axis=1)


def add_all_scores(adata):
    """
    Adds all scores listed above.
    :param adata: AnnData objects, contains genes listed in markers in adata.X
    :return: results in adata.obs[name] for all markers
    """
    for name, genes in marker_collection.items():
        add_score(adata, genes, name)


def add_tc_scores(adata):
    """
    Adds the scores from Szabo et al.
    https://www.nature.com/articles/s41467-019-12464-3
    :param adata: AnnData objects, contains genes listed in markers in adata.X
    :return: results in adata.obs[name] for all markers
    """
    tcell_signatures = pd.read_excel('../../data/scores/41467_2019_12464_MOESM7_ESM.xlsx', skiprows=1)
    for score in tcell_signatures.columns:
        add_score(adata, tcell_signatures[score].values, f'{score}_score')


def hash_solo_by_sample(adata, column_sample, samples, list_hash_cols, col_name, n_noise_barcodes):
    adata.obs[col_name] = 'NaN'

    dfs_donor = []
    for s, hashtag_cols in zip(samples, list_hash_cols):
        adata_tmp = adata[adata.obs[column_sample] == s].copy()
        adata_tmp.obs = adata_tmp.obs.drop(col_name, axis=1)

        sc.external.pp.hashsolo(adata_tmp, hashtag_cols, number_of_noise_barcodes=n_noise_barcodes)
        adata_tmp.obs = adata_tmp.obs.rename(columns={'Classification': col_name})

        dfs_donor.append(adata_tmp.obs[[col_name]])

    dfs_donor = pd.concat(dfs_donor)
    dfs_donor = dfs_donor.loc[adata.obs.index]
    adata.obs[col_name] = dfs_donor


def extract_clonotype_information(adata, information, name):
    clono_info = adata.obs[[f'IR_VJ_1_{information}', f'IR_VJ_2_{information}',
                            f'IR_VDJ_1_{information}', f'IR_VDJ_2_{information}',
                            'clone_id']].copy()
    for col in clono_info.columns:
        if col == 'clone_id':
            continue
        clono_info[col] = clono_info[col].astype(str)
        clono_info[col] = clono_info[col].replace('nan', np.nan)
    clono_info = clono_info.groupby('clone_id')

    def assign_clonotype_seqs(chains):
        if len(chains) == 1:
            alpha = chains[f'IR_VJ_1_{information}'].iloc[0]
            beta = chains[f'IR_VDJ_1_{information}'].iloc[0]
            return f'{alpha} {beta}'
        beta_1 = chains[f'IR_VDJ_1_{information}'].value_counts()
        beta_joint = pd.concat([chains[f'IR_VDJ_1_{information}'], chains[f'IR_VDJ_2_{information}']]).value_counts()

        if len(beta_joint) > 1 and beta_joint.values[0] == beta_joint.values[1]:
            beta = beta_1.keys()[0]
        else:
            beta = beta_joint.keys()[0]

        alpha_1 = chains[f'IR_VJ_1_{information}'].value_counts()
        alpha_joint = pd.concat([chains[f'IR_VJ_1_{information}'], chains[f'IR_VJ_2_{information}']]).value_counts()

        if len(alpha_joint) > 1 and alpha_joint.values[0] == alpha_joint.values[1]:
            alpha = alpha_1.keys()[0]
        else:
            alpha = alpha_joint.keys()[0]
        return f'{alpha} {beta}'

    clone_dict = {ct: assign_clonotype_seqs(chains) for ct, chains in clono_info}
    clone_dict['nan'] = 'nan'
    adata.obs[name] = adata.obs['clone_id'].map(clone_dict)
    adata.uns[f'{name}_dict'] = clone_dict


def extract_nucleotide_seq(adata, col_aa):
    def get_majority_nt(df, chain_id):
        aa_seq = df.name.split(' ')
        if chain_id == 'VJ':
            aa_seq = aa_seq[0]
        else:
            aa_seq = aa_seq[1]
        chain_1 = df[f'IR_{chain_id}_1_junction'].values[df[f'IR_{chain_id}_1_junction_aa'].values == aa_seq]
        chain_2 = df[f'IR_{chain_id}_2_junction'].values[df[f'IR_{chain_id}_2_junction_aa'].values == aa_seq]
        all_chains = np.concatenate([chain_1, chain_2])
        return Counter(all_chains).most_common(1)[0][0]
    df_base = adata[adata.obs['clone_id'] != 'nan'].obs

    aa_2_nt_alpha = df_base.groupby(col_aa).progress_apply(lambda x: get_majority_nt(x, 'VJ'))
    aa_2_nt_beta = df_base.groupby(col_aa).progress_apply(lambda x: get_majority_nt(x, 'VDJ'))
    aa_2_nt = aa_2_nt_alpha + ' ' + aa_2_nt_beta
    aa_2_nt = dict(aa_2_nt)
    aa_2_nt['nan'] = 'nan'
    adata.obs[f'{col_aa}_nt'] = adata.obs[col_aa].map(aa_2_nt)


def assign_mait(adata, gene_combination=True, cdr3=True):
    def mait_by_row(row):
        if gene_combination:
            if 'TRAJ33' in row['j_genes'] and 'TRAV1-2' in row['v_genes']:
                if 'TRBV20-1' in row['v_genes'] or 'TRBV6' in row['v_genes']:
                    return 'True'
        if cdr3:
            if 'CAVMDSSYKLIF' in row['clonotype_sequence']:
                return 'True'
        return 'False'
    adata.obs['has_MAIT'] = adata.obs.apply(mait_by_row, axis=1)
