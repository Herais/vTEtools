import numpy as np
import pandas as pd
import scipy as sp
import scanpy as sc
from anndata import AnnData

class SCRNA(object):
    """
    """
    def __init__(self):
        """
        Parameters
        ----------

        """
        super(SCRNA, self).__init__()

    
    @staticmethod
    def process_df_to_scanpy(df):

        # filter columns for TE(s) and marker(s)    
        data = df.copy()
        features = data.columns
        cells = data.index
        
        import scipy as sp
        data = sp.sparse.csr_matrix(data.to_numpy())
        data.astype('float32')
        
        ad = AnnData(data, obs={'obs_names': cells}, var={'var_names': features})
        del data
        return ad

    @staticmethod
    def read_pl_rank_genes_groups_to_df(
            adata,
            cols=['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']):

        ls = []
        for col in cols:
            ls.append(adata.uns['rank_genes_groups'][col])

        df = pd.DataFrame(ls).T
        df = df.map(lambda x: x[0])
        df.columns = cols

    return df.copy()

    @staticmethod
    def screenTE_adata2(adata,
                        pvals=None,
                        pvals_adj=None,
                        min_logfoldchange=None,
                        method='wilcoxon',
                        layer='norm_to_TE_include'):

        #pair = (reference, group)
        pairs = [('Sham', 'TAC'),
                (['TAC', 'TACJQ1']),
                (['TACJQ1', 'TACJQ1w']),
                ]

        xy2df = {}

        if pvals:
            _pvals = pvals
        else:
            _pvals = 1

        if pvals_adj:
            _pvals_adj = pvals_adj
        else:
            _pvals_adj = 1

        if min_logfoldchange:
            _min_logfoldchange = abs(min_logfoldchange)
        else:
            _min_logfoldchange = 0

        if method:
            _method = method
        else:
            _method = 't-test'

        for pair in pairs:
            xy = "{} vs. {}".format(pair[1], pair[0])
            xy2df[xy] = {}
            if _method == 'wilcoxon':
                sc.tl.rank_genes_groups(adata,
                                        groupby='condition',
                                        groups=[pair[1]],
                                        reference=pair[0],
                                        method=_method,
                                        tie_correct=True,
                                        layer=layer)
            elif _method == 't-test':
                sc.tl.rank_genes_groups(adata,
                                        groupby='condition',
                                        groups=[pair[1]],
                                        reference=pair[0],
                                        method=_method,
                                        layer=layer)

            xy2df[xy]['df'] = read_pl_rank_genes_groups_to_df(adata)

            xy2df[xy]['up'] = xy2df[xy]['df'][xy2df[xy]['df'].logfoldchanges > 0]
            xy2df[xy]['up'] = xy2df[xy]['up'][xy2df[xy]['up'].pvals <= _pvals]
            xy2df[xy]['up'] = xy2df[xy]['up'][xy2df[xy]
                                            ['up'].pvals_adj <= _pvals_adj]
            xy2df[xy]['up'] = xy2df[xy]['up'][xy2df[xy]
                                            ['up'].logfoldchanges > _min_logfoldchange]

            xy2df[xy]['down'] = xy2df[xy]['df'][xy2df[xy]['df'].logfoldchanges < 0]
            xy2df[xy]['down'] = xy2df[xy]['down'][xy2df[xy]['down'].pvals <= _pvals]
            xy2df[xy]['down'] = xy2df[xy]['down'][xy2df[xy]
                                                ['down'].pvals_adj <= _pvals_adj]
            xy2df[xy]['down'] = xy2df[xy]['down'][xy2df[xy]
                                                ['down'].logfoldchanges < -_min_logfoldchange]

        return xy2df
