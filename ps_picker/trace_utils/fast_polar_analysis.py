# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# fast_polar_analysis.m

import numpy as np
import scipy.linalg as la


def fast_polar_analysis(times, compute_wind, tracex, tracey, tracez,
                        analyze_wind=4.):
    """
    Compute the polarity parameters of a 3 component seismogram.

    Only computes around user specified times since this analysis is time
    consuming.  Other values are set to zero
    :param times: list of times around which to compute polarity
    :param compute_wind: seconds around times to analyze
    :param trace?: traces
    :param analyze_wind = length of analysis window for each calculation
        (seconds)
    :returns: rect, azi, dip
    """
    # fast_polar_analysis.m:17
    sr = tracex.stats.sampling_rate
    assert tracey.stats.sampling_rate == sr,\
        "tracey has a different sampling rate from tracex"
    assert tracez.stats.sampling_rate == sr,\
        "tracez has a different sampling rate from tracex"

    # Create a list of all indices to investigate: the specified compute
    # window around all picks minus any duplicate samples or any that are
    # too close to the window edges for a proper analysis
    n_half_compute = round(sr * compute_wind / 2)
    nsamps_min = np.min([len(tracex.data), len(tracey.data), len(tracez.data)])
    ind_vec = []
    # for i in arange(1, length(indices)).reshape(-1):
    for t in times:
        index = tracex.times('utcdatetime').searchsorted(t)
        iwind = np.arange(index - n_half_compute, index + n_half_compute)
        ind_vec = np.concat([ind_vec, iwind])

    n_half_analyze = round(analyze_wind * sr / 2)
    # Get rid of duplicate or out-of-bounds values
    # ind_vec[ind_vec < 0] = []
    # ind_vec[ind_vec >= nsamps_min] = []
    ind_vec[ind_vec - n_half_analyze <= 0] = []
    ind_vec[ind_vec + n_half_analyze > nsamps_min] = []
    ind_vec = sorted(ind_vec)
    ind_vec[np.diff(ind_vec) == 0] = []

    rectilinP = tracex.copy()
    rectilinP.data = np.zeros(nsamps_min)
    azimuthP = rectilinP.copy()
    dipP = rectilinP.copy()

    # define start and end time for the analysis
    # loop over moving time window
    # for k in concat([ind_vec.T]).reshape(-1):
    for k in ind_vec:
        xP = tracex.data[k - n_half_analyze: k + n_half_analyze]
        yP = tracey.data[k - n_half_analyze: k + n_half_analyze]
        zP = tracez.data[k - n_half_analyze: k + n_half_analyze]
        MP = np.cov(np.array([xP, yP, zP]))
        D, V = la.eig(MP)
        i_sort = np.argsort(D)  # sort from min (0) to max (2)
        D = D[i_sort]
        V = V[:, i_sort]
        # 2 = major axis, 0 = minor axis, 1 = intermediate
        rectilinP.data[k] = 1 - ((D[0] + D[1]) / (2 * D[2]))
        V_major = V[:, 2]
        azimuthP.data[k] = np.degrees(np.atan(V_major[1] / V_major[0]))
        dipP.data[k] = np.degrees(np.atan(V_major[2]
                                  / np.sqrt(V_major[1]**2 + V_major[0]**2)))

#         rectilinP=[NaN(n_half_analyze,h); rectilinP];
#         azimuthP=[NaN(n_half_analyze,h); azimuthP];
#         dipP=[NaN(n_half_analyze,h); dipP];
    return rectilinP, azimuthP, dipP


if __name__ == '__main__':
    pass
