import numpy as np
from scipy.signal import lfilter
from ..logger import log

from obspy.core.stream import Stream as obspy_Stream
from obspy.core.stream import Trace


def smooth_filter(traces, n_smooth):
    """
    Smooth a trace or a list of traces

    Uses the scipy.signal.lfilter() function but, if the trace begins with NaN
    values, it smooths the entry into "good" data and does not erase the
    first non-nan values.

    Reproduces Christian's code, which always ramps up the beginning,
    even if there are no NaNs

    :param traces_in: trace of list of traces
    :param n_smooth: size of the smoothing window in SAMPLES
    :returns: smoothed data
    """
    n_smooth = int(n_smooth)
    bare_trace = False
    of isinstance(traces, Trace):
        traces = [traces]
        bare_trace = True
    assert isinstance(traces, list) or isinstance(traces, obspy_Stream):
    # Eliminate traces filled with nan
    cleaned_traces_in = [t for t in traces if not np.all(np.isnan(t.data))]

    smoothed_traces = []
    for tr in cleaned_traces_in:
        i_nan = np.nonzero(np.isnan(tr.data))[0]
        if len(i_nan) > 0:
            tr.data[i_nan] = 0
        # tr.data[np.isnan(tr.data)] = 0
        smooth_tr = tr.copy()
        smooth_tr.data = lfilter(np.divide(np.ones(n_smooth), n_smooth),
                                 1., tr.data)
        correction_vec = np.ones(len(tr.data))
        # Find first non-nan value in tr and ramp up over n_smooth
        ind_correc = np.nonzero(np.isfinite(tr.data))[0]
        if len(ind_correc) > 0:
            correction_vec[ind_correc[0]: ind_correc[0] + n_smooth] =\
                np.divide(np.arange(0, n_smooth), n_smooth)
        smooth_tr.data *= correction_vec
        i_nan = np.nonzero(np.isnan(tr.data))[0]
        if len(i_nan) > 0:
            smooth_tr.data[i_nan] = np.nan
        smoothed_traces.append(smooth_tr)
    if len(smoothed_traces) == 0:
        return None
    if bare_trace:
        return smoothed_traces[0]
    else:
        return smoothed_traces


if __name__ == '__main__':
    pass
