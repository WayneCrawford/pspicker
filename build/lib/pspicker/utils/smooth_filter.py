import numpy as np
from scipy.signal import lfilter
# from ..logger import log

from obspy.core.stream import Stream
from obspy.core.stream import Trace


def smooth_filter(traces, n_smooth):
    """
    Smooth a trace or a list of traces

    Uses the scipy.signal.lfilter() function

    :param traces_in: trace of list of traces
    :param n_smooth: size of the smoothing window in SAMPLES
    :returns: smoothed data
    """
    n_smooth = int(n_smooth)
    bare_trace = False
    if isinstance(traces, Trace):
        traces = [traces]
        bare_trace = True
    assert isinstance(traces, list) or isinstance(traces, Stream)

    smoothed_traces = []
    for tr in traces:
        smooth_tr = tr.copy()
        smooth_tr.data = lfilter(np.divide(np.ones(n_smooth), n_smooth),
                                 1., tr.data)
        smoothed_traces.append(smooth_tr)
    if len(smoothed_traces) == 0:
        return None
    if bare_trace:
        return smoothed_traces[0]
    return smoothed_traces


if __name__ == '__main__':
    pass
