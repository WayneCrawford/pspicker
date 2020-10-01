from obspy.core.trace import Trace


def mean_trace(traces):
    """
    Calculate the mean trace

    :param traces: stream or list of traces
    :returns: trace object
    """
    if isinstance(traces, Trace):
        return traces
    assert isinstance(traces[0], Trace), "traces[0] is not a trace"
    if len(traces) == 1:
        return traces[0].copy

    data_len = len(traces[0].data)
    for tr in traces[1:]:
        assert len(tr.data) == data_len, 'traces are not the same length'

    mean_trace = traces[0].copy
    for tr in traces[1:]:
        mean_trace.data += tr.data
    mean_trace.data /= len(traces)
    return mean_trace


if __name__ == '__main__':
    pass
