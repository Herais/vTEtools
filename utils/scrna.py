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