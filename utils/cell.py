import pandas as pd
import numpy as np

class CELL(object):
    """
    """
    def __init__(self):
        """
        Parameters
        ----------

        """
        super(CELL, self).__init__()

    def get_celltype(df, celltype2markers):
        """
        """
        markers = list(itertools.chain.from_iterable(celltype2markers.values()))
        markers = list(set(markers))
        df = df.loc[:,df.columns.isin(markers)].copy()
        cols = []
        for k, v in celltype2markers.items():
            df['ct_{}'.format(k)] = df.loc[:, df.columns.isin(v)].sum(axis=1)
            cols.append('ct_{}'.format(k))
        df['celltype'] = df[cols] \
            .apply(lambda x: cols[x.argmax()], axis=1) \
            .apply(lambda x: x[3:])
        df['cellassignment'] = df[cols].sum(axis=1) > 0
        df.loc[~df.cellassignment, 'celltype'] = 'other'
        
        return df['celltype'].copy(), df.copy()
    
    def assign_celltype_strict(df, celltype2markers, relax=0):
    
        markers = list(itertools.chain.from_iterable(celltype2markers.values()))
        markers = list(set(markers))
        df = df.loc[:,df.columns.isin(markers)].copy()
        
        for celltype, markers in celltype2markers.items():
            num_unique_markers = len(markers)
            df['celltype'] = 'other'
            df_celltype = df.loc[:,df.columns.isin(markers)]
        
            S_num_unique_markers_identified = df_celltype[df_celltype > 0].sum(axis=1)
            idx_all_markers_present = df_celltype \
                .loc[(S_num_unique_markers_identified == num_unique_markers-relax), :] \
                .index
            df['{}_w_allmarkers'.format(celltype)] = False
            df.loc[idx_all_markers_present, '{}_w_allmarkers'.format(celltype)] = True

            df['{}_markers_count'.format(celltype)] = df_celltype.sum(axis=1)
        
        S_cells_num_complete_markers = df \
            .loc[:, df.columns.str.contains('_w_allmarkers')].sum(axis=1)

        df_cells_w_complete_markers = df \
            .loc[S_cells_num_complete_markers > 0, df.columns.str.contains('_w_allmarkers')]
        df_cells_w_complete_markers.columns = \
            [col.split('_')[0] for col in df_cells_w_complete_markers.columns]
        
        df_cells_w_complete_markers_count = df \
            .loc[S_cells_num_complete_markers > 0, df.columns.str.contains('_markers_count')]
        df_cells_w_complete_markers_count.columns = \
            [col.split('_')[0] for col in df_cells_w_complete_markers_count.columns]
        
        df_cells_w_complete_markers_count_masked = \
            df_cells_w_complete_markers_count \
                .mask(~df_cells_w_complete_markers).fillna(0)
        
        cols = df_cells_w_complete_markers_count_masked.columns
        S_cells_with_celltype = df_cells_w_complete_markers_count_masked \
            .apply(lambda x: cols[x.argmax()], axis=1)
        df.loc[S_cells_with_celltype.index, 'celltype'] = S_cells_with_celltype

        return df.loc[:, 'celltype'].copy()