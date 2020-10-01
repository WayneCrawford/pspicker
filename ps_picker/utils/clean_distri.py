import numpy as np


def clean_distri(values, n_std=3, mode='median', min_vals=3):
    """
    Remove extreme values in a vector

    :param values: input vector
    :param n_std: erase values for which x is more than n_std standard
        deviations from mode(x)
    :param mode: 'mean' or 'median' [ default = 'median' ]
    :param min_vals:  do not clean if vector is shorter than this[default=3]
    :returns: list of values after cleaning, indices of values retained

    >>> clean_distri([0, 3, 4, 6, 3, 5, 25])
    (array([0, 3, 4, 6, 3, 5]), array([0, 1, 2, 3, 4, 5]))
    """
    if len(values) < min_vals:
        return values, list(range(len(values)))

    x = np.array(values)
    if mode == 'mean':
        dev = x - np.mean(x)
    else:
        dev = x - np.median(x)

    # Calculate std of values
    lim = round(0.68*len(x))  # assume gaussian around median
    if lim == 0:
        lim = 1
    std = np.sort(np.abs(dev))[lim]

    in_rm = np.nonzero(np.abs(dev) <= n_std * std)[0]

    return x[in_rm], in_rm


if __name__ == "__main__":
    import doctest
    doctest.testmod()
