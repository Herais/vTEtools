import numpy as np
import pandas as pd
import matplotlib
from matplotlib import cm
import seaborn as sns
import matplotlib.pyplot as plt

def generate_volcano_plot(df,
                          colx: str,
                          coly: str,
                          title: str = None,
                          vline_left_FC=-0.5,
                          vline_right_FC=0.5,
                          p_threshold=0.05,
                          cmap='coolwarm',
                          figsize=(4, 4),
                          ):

    ret = {}

    # assign colors
    C = df[colx]
    C = C.apply(lambda x: -1 if x < 0 else x)
    C = C.apply(lambda x: 0 if x > 0 else x)
    C = C.apply(lambda x: 0 if x == 0 else x)

    hline_p_threshold = df.loc[df.pvals < p_threshold, 'neglog10pvals'].min()

    fig, ax = plt.subplots()
    df.plot(kind='scatter',
            x=colx,
            y=coly,
            title=title,
            s=1,
            c=C,
            cmap=cmap,
            ax=ax,
            figsize=figsize,
            )
    plt.vlines(vline_left_FC, ymin=df[coly].min(
    ), ymax=df[coly].max(), color='blue', linestyles='dashed')
    plt.vlines(vline_right_FC, ymin=df[coly].min(
    ), ymax=df[coly].max(), color='red', linestyles='dashed')
    plt.vlines(0, ymin=df[coly].min(), ymax=df[coly].max(),
               color='black', linestyles='dotted')
    plt.hlines(hline_p_threshold, xmin=df[colx].min(), xmax=df[colx].max(),
               color='black', linestyles='dashed', label='p={}'.format(p_threshold))

    plt.legend()
    cbar = ax.collections[0].colorbar
    cbar.remove()

    ret['fig'] = fig
    ret['ax'] = ax
    ret['upregulated'] = df[(df.pvals<= p_threshold)
                            & (df.logfoldchanges > vline_right_FC)]
    ret['downregulated'] = df[(df.pvals <= p_threshold)
                              & (df.logfoldchanges < vline_left_FC)]

    return ret
