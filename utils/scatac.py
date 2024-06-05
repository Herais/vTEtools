import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp
from anndata import AnnData

class SCATAC(object):
    """
    """
    def __init__(self):
        """
        Parameters
        ----------

        """
        super(SCATAC, self).__init__()
    
    @staticmethod
    def process_df_to_scanpy(df):

        # filter columns for TE(s) and marker(s)    
        data = df.copy()
        features = data.columns
        cells = data.index
            
        data = sp.sparse.csr_matrix(data.to_numpy())
        data.astype('float32')

        '''
        oh = open('gene_names.{0}.tsv'.format(os.path.split(filename)[1]), 'w')
        for g in genes:
            oh.write('%s\n' % g)
        oh.close()
        '''
        
        print('Loaded {0}'.format(name))
        ad = AnnData(data, obs={'obs_names': cells}, var={'var_names': features})
        del data
        return ad
    
    @staticmethod
    def read_pl_rank_genes_groups_to_df(adata):

        cols = ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']
        uns_col = 'rank_genes_groups'

        ls = []
        for col in cols:
            if col in adata.uns[uns_col].keys():
                ls.append(adata.uns[uns_col][col])

        df = pd.DataFrame(ls).T
        df = df.map(lambda x: x[0])
        df.columns = cols

        df = pd.concat([df.set_index('names'),
                        adata.uns[uns_col]['pts'].add_prefix('pts_')],
                       axis=1)

        print(adata.uns[uns_col]['params'])

        return df.copy()

    @staticmethod
    def read_pl_rank_genes_groups_filtered_to_df(adata, uns_col='rank_genes_groups_filtered'):

        cols = ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']

        ls = []
        for col in cols:
            if col in adata.uns[uns_col].keys():
                ls.append(adata.uns[uns_col][col])

        df = pd.DataFrame(ls).T
        df = df.map(lambda x: x[0])
        df.columns = cols

        df = df.dropna(subset=['names']).set_index('names')
        df_pts = adata.uns[uns_col]['pts'].add_prefix('pts_')
        df[['pts_TAC', 'pts_Sham']] = df_pts

        print(adata.uns[uns_col]['params'])

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

            mySCRNA = SCRNA()

            xy2df[xy]['df'] = mySCRNA.read_pl_rank_genes_groups_to_df(adata)

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

    @staticmethod
    def generate_pattern_descriptions(xy2df, pairchanges, patterns):

        pattern2df = {}
        
        for pattern in patterns:
            
            pattern_name = '-'.join(pattern)
            pattern2df[pattern_name] = {}
            
            lsls = []
            for i, move in enumerate(pattern):
                xy = pairchanges[i]
                pair_pattern = ("{} {}".format(xy, move.upper()))
        
                
                df= xy2df[xy][move]
                pattern2df[pattern_name][pair_pattern] = df['names'].shape[0]
                
                lsls.append(df['names'].to_list())
        
            common = set(lsls[0])
            for ls in lsls[1:]:
                common = common.intersection(ls)
            pattern2df[pattern_name]['common'] = list(common)
            pattern2df[pattern_name]['num_common'] = len(common)
        return pattern2df