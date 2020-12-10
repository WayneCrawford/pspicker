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
        cluster_data = fclusterdata(np.array([[p.time for p in picks]]).T,
                                    t=window_sec,
                                    criterion='distance')
        cluster_groups = dict()
        for cnum in cluster_data:
            cluster_groups[cnum] = cluster_groups.get(cnum, 0) + 1
        max_length = max(cluster_groups.values())
        max_keys = [k for k, v in cluster_groups.items() if v==max_length]
             
        if len(max_keys) > 1:
            warnings.warn('{} clusters of max_length ({:d}), returning first'
                          .format(len(max_keys), max_length))
        picks = [p for p, i in zip(picks, cluster_data) if i==max_keys[0]]
    return picks


if __name__ == '__main__':
    pass
