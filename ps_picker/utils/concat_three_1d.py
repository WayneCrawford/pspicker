import numpy as np


def concat_three_1d(left, center, right):
    """
    Concatenate three 1d ndarrays, of which the left and right may be empty
    """
    # print(len(left), len(center), len(right))
    return np.concatenate((left, center, right))
    # if len(left) == 0:
    #     if len(right) == 0:
    #         return center
    #     else:
    #         return np.concatenate((center, right))
    # else:
    #     if len(right) == 0:
    #         return np.concatenate((left, center))
    #     else:
    #         return np.concatenate((left, center, right))
