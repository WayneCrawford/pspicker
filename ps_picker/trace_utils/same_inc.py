import numpy as np
from obspy.core.trace import Trace


def same_inc(traces_old, firsttime, lasttime):
    """
    Return a section of an input trace and NaNs elsewhere

    :param traces_old: trace or list of input traces
    :param firsttime: time of first sample to keep
    :param lasttime: time of last sample to keep
    :returns: output trace or list of traces
    """
    # same_inc.m:6
    bare_trace = False
    if isinstance(traces_old, Trace):
        traces_old = [traces_old]
        bare_trace = True
    if len(traces_old) == 0:
        return []

    traces_new = []
    for trace in traces_old:
        tr = trace.copy()
        sr = tr.stats.sampling_rate
        stime = tr.stats.starttime
        start_ind = int(np.floor((firsttime - stime) * sr))
        end_ind = int(np.floor((lasttime - stime) * sr))
        if start_ind < 0:
            start_ind = 0
        if end_ind >= len(tr.data):
            end_ind = len(tr.data) - 1
        if start_ind > end_ind:
            start_ind = end_ind
        length = len(tr.data)
        left_mat = np.empty(start_ind)
        right_mat = np.empty(length - (end_ind + 1))
        left_mat[:] = np.nan
        right_mat[:] = np.nan
        tr.data = np.concatenate((left_mat, tr.data[start_ind: end_ind],
                                  right_mat))
        traces_new.append(tr)

    if bare_trace:
        return traces_new[0]
    else:
        return traces_new


if __name__ == '__main__':
    pass
