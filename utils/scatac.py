import numpy as np
import pandas as pd
import scanpy as sc
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