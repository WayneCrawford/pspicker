# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# cluster_cleaning.m
import numpy as np
from scipy.cluster.hierarchy import fclusterdata


def cluster_clean(window_sec=None, picks=None):
    """
    Return indices of picks that do not fit in the cluster

    :param window_sec: cluster window length
    :param picks: input picks

    :returns: picks that fit the cluster criteria
    """

    if len(picks) > 1:
        clust_group = fclusterdata(np.array([[p.time for p in picks]]).T,
                                   t=window_sec,
                                   criterion='distance')
        print(clust_group)
        # for kk in arange(1,max(clust_group)).reshape(-1):
        lengths = [len(x) for x in clust_group]
        i_clust_max = np.argmax(lengths)
        if len(i_clust_max) > 1:
            picks = []  # Don't even chose one?
        else:
            picks = picks[clust_group[i_clust_max]]
    return picks


if __name__ == '__main__':
    pass
