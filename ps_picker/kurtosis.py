"""
Kurtosis-specific routines
"""
import warnings

import numpy as np
from scipy.signal import lfilter
from obspy.core.stream import Stream
from obspy.core.trace import Trace

from .trace_utils import smooth_filter


class Kurtosis():
    """
    Class for picking seismograms using the Kurtosis
    """
    def __init__(self, freq_bands, wind_lengths, n_smooth):
        """
        :param FBs: list of [min, max] frequency bands
        :param window_lengths: list of window lengths (seconds)
        :param n_smooth: smoothings to apply in calculating Kurtosis (samples)
        """
        assert isinstance(freq_bands, list), "freq_bands is not a list"
        assert isinstance(wind_lengths, list), "wind_lengths is not a list"
        assert isinstance(n_smooth, int), "n_smooth is not an int"

        self.freq_bands = freq_bands
        self.wind_lengths = wind_lengths
        self.n_smooth = n_smooth
        self.extrem_smooths = None  # smoothings for extrema following
        
        # Mean traces over freq_bands, wind_lengths & smoothing
        self.mean_kurtosis = None
        self.mean_cumulative_kurtosis = None
        
        # Gradients of the mean_cumulative kurtosis, for each extrem_smoothing
        self.kurto_gradients = None

    def pick_trace(self, trace, n_picks, starttime=None, endtime=None,
                   extrem_type='mini', extrem_smooths=[1, 5, 10],
                   extrem_normalize=False, extrem_which=None):
        """
        Pick a trace using the Kurtosis

        :param trace: waveform trace
        :param n_picks: number of picks to return
        :param starttime: first time of interest
        :param endtime: last time of interest
        :param extrem_type: 'mini' or 'maxi', depending on wheter you want to
            follow minima or maxima
        :param extrem_smooths: smoothing to apply when following extrema
        :param extrem_normalize: normalize the cumulative Kurtosis gradient
        :param which:
            'first': select first 'n_follow' extrema > 0.1, order right to left
            otherwise: select the 'n_follow' biggest extrema, order biggest
                       to smallest
        :returns:  indices of extrema,
                   values of extrema
        """
        self.extrem_smooths = extrem_smooths
        self.trace2kurtocum(trace, starttime, endtime)
        t, v = self.follow_extrem(type_=extrem_type, n_follow=n_picks,
                                  normalize=extrem_normalize,
                                  sense=extrem_which)
        return t, v

    def trace2kurtocum(self, trace, starttime=None, endtime=None, debug=False):
        """
        Calculate cumulative kurtosis over bandwidths and windows

        FW stands for Frequency Window because it computes for
        several frequency bandwidth and window sizes

        Puts mean kurtosis over all windows and freq_bands in
        self.mean_cumulative_kurtosis
        :param trace: the raw trace (works one trace, not a stream)
        :param FBs: list of [min, max] frequency bands
        :param window_lengths: list of window lengths (seconds)
        :param starttime: is the first time of interest
        :param endtime: is the last time of interest
        """
        # trace2FWkurto.m:26
        assert isinstance(trace, Trace), "trace is not an obspy Trace"

        sr = trace.stats.sampling_rate
        data_length = trace.stats.endtime - trace.stats.starttime
        if starttime is None:
            starttime = trace.stats.starttime
        elif starttime < trace.stats.starttime:
            starttime < trace.stats.starttime
        if endtime is None:
            endtime = trace.stats.endtime
        elif endtime > trace.stats.endtime:
            endtime = trace.stats.endtime

        # Filter traces in different frequency bands
        B = []
        for FB in self.freq_bands:
            f = trace.copy()
            if debug:
                print(f'trace2FWkurto: filtering from {FB[0]} to {FB[1]} Hz')
            f.detrend('demean')
            f.filter(type='bandpass', freqmin=FB[0], freqmax=FB[1], corners=3)
            f = f.slice(starttime, endtime)
            f.data[np.arange(0, 50)] = 0
            B.append(f)

        # 2-level lists: : 1st dim: window_lengths, 2nd dim: freq bands
        K = []  # all kurtosises:
        C = []  # all cumulative, detrended kurtosises
        for win_len in self.wind_lengths:
            if win_len > data_length:
                warnings.warn(
                    'Kurtosis window longer than data window '
                    '({:2g}s > {:4g}s), skipping!'.format(win_len,
                                                          data_length / sr))
            else:
                win_samps = int(np.floor(win_len * sr)) + 1
                if win_samps > len(trace.data):
                    win_samps = len(trace.data)
                kurtos = []
                for t in B:
                    k = _fast_kurtosis(t, win_samps)
                    if debug:
                        print('plotting trace, kurtosis for win_samps='
                              f'{win_samps}, f = {self.freq_bands[0]}')
                        k.stats.channel = 'KUR'
                        Stream([k, t]).plot(equal_scale=False)
                    kurtos.append(k)
                f_kurto_cums = smooth_filter(kurtos, self.n_smooth)
                kurto_cums, _ = _f_cumul(f_kurto_cums)
                lines = _f_segment(kurto_cums)
                corr_kurto_cums = kurto_cums.copy()
                if debug:
                    print(f'trace2FWkurto: plotting kurt_cums & fit lines')
                for c, l in zip(corr_kurto_cums, lines):
                    if debug:
                        l.stats.channel = 'LIN'
                        Stream([c, l]).plot()
                    c.data -= l.data
                C.append(corr_kurto_cums)
                K.append(kurtos)
        self.mean_kurtosis = _mean_trace(K)
        self.mean_cumulative_kurtosis = _mean_trace(C)

    def follow_extrem(self, type_='mini', n_follow=2,
                      normalize=False, sense=None, debug=False):
        """
        Find extremas using smoothed versions of the mean cumulative kurtosis

        Steps:
          1) locate the 'n' first extremas of the smoothest function
          2) refine these points step by step through less and less smooth
             functions

        Uses self.mean_cumulative kurtosis
        :param type: 'mini' or 'maxi':  follow minima or maxima
        :param n_follow: number of extrema to follow
        :param normalize: normalize the gradient?
        :param sense:
            'first': select first 'n_follow' extrema > 0.1, order right to left
            otherwise: select the 'n_follow' biggest extrema, order biggest
                       to smallest
        :returns:  UTCDateTimes of extrema,
                   values of extrema
        """
        # follow_extrem.m:29
        # Parameters
        assert not len(self.mean_cumulative_kurtosis) == 0,\
            'no mean cumulative kurtosis!'
        data = self.mean_cumulative_kurtosis
        st = self.mean_cumulative_kurtosis.stats.starttime
        et = self.mean_cumulative_kurtosis.stats.endtime
        
        gradients = []
        # Put the strongest smoothing first
        self.extrem_smooths.sort(reverse=True)
        # for i in arange(1,length(extrem_smooths)).reshape(-1):
        for smoothing in self.extrem_smooths:
            v_smooth = smooth_filter(data, smoothing)
            v_smooth = _cum2grad(v_smooth, normalize)
            ext_indices, _ = _loca_ext(v_smooth, st, et, type_)
            if debug:
                print(f'smoothing={smoothing}: i_extrema={ext_indices}')
            # v_smooth.plot()
            gradients.append({'indices': ext_indices, 'trace': v_smooth})

        # create a list of {index:  value:} dictionaries with all of the
        # detected extrema in the smoothest trace
        all_extrema = [{'index': i, 'value': gradients[0]['trace'].data[i]}
                       for i in gradients[0]['indices']]
        # print(f'all_extrema = {all_extrema}')
        # print(f'gradients[0] = {gradients[0]}')
        # Choose up to 'n_follow' extrema in the smoothest trace
        if sense == 'first':
            # If extrema are tiny, take the biggest
            big_extrema = [x for x in all_extrema if x['value'] >= 0.1]
            if len(big_extrema) > 1:
                ext = sorted(big_extrema, key=lambda k: k['index'])
                try:
                    b = [x for x in ext[:n_follow:-1]]
                except Exception:
                    b = [x for x in ext[::-1]]
            else:
                ext = sorted(all_extrema, key=lambda k: np.abs(k['value']),
                             reverse=True)
                try:
                    b = [x for x in ext[:n_follow]]
                except Exception:
                    b = [x for x in ext]
        else:
            ext = sorted(all_extrema, key=lambda k: np.abs(k['value']),
                         reverse=True)
            try:
                b = [x for x in ext[:n_follow]]
            except Exception:
                b = [x for x in ext]
        if len(b) == 0:
            return [x['trace'] for x in gradients], b, []

        # Sharpen the indices/values using the smaller smoothing values
        # for k in arange(2, n).reshape(-1):
        for M, smoothing in zip(gradients[1:], self.extrem_smooths[1:]):
            if len(M['indices']) == 0:
                warnings.warn(
                    f'No extrema found for {smoothing:d}-sample '
                    'smoothing of Kurtosis')
            else:
                # for j in arange(1, len(b)).reshape(-1):
                c = []
                # print(f'b = {b}')
                for extrem in b:
                    # print(f'extrem = {extrem}')
                    # print(f'M = {M}')
                    differ = abs(extrem['index'] - M['indices'])
                    if np.any(differ <= 40):
                        i = np.argmin(differ)
                        c.append({'index': M['indices'][i],
                                  'value': M['trace'].data[i]})
                    else:
                        c.append(extrem)
                b = c.copy()
        self.kurto_gradients = [x['trace'] for x in gradients]

        sr = self.mean_cumulative_kurtosis.stats.sampling_rate
        return [st + x['index']/sr for x in b], [x['value'] for x in b]


def _fast_kurtosis(trace, win_samps):
    """
    Compute kurtosis really quickly using "filter" function

    could I just use obspy kurtosis?

    :param trace: one trace
    :param win_samps: number of samples in the sliding window
    :returns: Kurtosis trace
    """
    assert isinstance(trace, Trace), "trace is not an obspy Trace"
    win_samps = int(round(win_samps))
    # fast_kurtosis.m:11
    if win_samps == 1:
        win_samps = 2

    # Set NaNs to 0 for computational stability
    f = trace.copy()
    f.detrend(type='demean')
    # f.data[np.isnan(f.data)] = 0

    # Compute kurtosis
    a = np.divide(np.ones(win_samps), float(win_samps))
    b = 1.
    m_2 = lfilter(a, b, f.data**2)
    m_4 = lfilter(a, b, f.data**4)
    out = trace.copy()
    out.data = np.divide(m_4, (m_2 ** 2))
    # Protect against edge effect
    out.data[:win_samps] = 0

    # Set any kurtosis value to nan for any indices within win_samples of
    # an NaN in the original data.  Not sure this does exactly like the
    # original
    for i in np.nonzero(trace.data == np.nan)[0]:
        if i:
            out.data[i: i + win_samps] = np.nan
    return out


def _f_cumul(f):
    """
    Calculate the positive gradient cumulative of f

    If the gradient is positive we take the cumulative, if the gradient is
    negative then, as long as the gradient is negative, the value is that
    of the last positive gradient. The output has then only positive
    gradients
         ___/
        /
    ___/
    :param f: stream or list of traces containing the data
    :returns: list of cumulative output traces, first value = 0
              list of cumulative output traces
    """
    # f_cumul.m:14
    bare_trace = False
    if isinstance(f, Trace):
        f = [f]
        bare_trace = True
    g = f.copy()
    for t in g:
        tdata = t.data.copy()    # Backup copy for info
        t.data[np.isnan(t.data)] = 0
        t.differentiate(method='gradient')
        t.data[t.data < 0] = 0
        t.data = np.cumsum(t.data)
        t.data[np.isnan(tdata)] = np.nan

    p = g.copy()

    # Subtract first non-nan value from g
    for t in g:
        ind = np.nonzero(np.isfinite(t.data))[0]
        if len(ind) > 0:
            t.data -= t.data[ind[0]]

    if bare_trace:
        return g[0], p[0]
    else:
        return g, p


def _f_segment(f):
    """
    Return a line segment between the first and last values of function.
    SHOULD BE REPLACEABLE BY DETREND('linear')

    Goes from first non-nan value to last non-nan value
    Input and output have the same size.
    :param f: stream or list of data traces
    :returns: list of line segment traces
    """
    bare_trace = False
    if isinstance(f, Trace):
        f = [f]
        bare_trace = True
    assert isinstance(f[0], Trace), f'f is list of {type(f[0])}s!'
    # f_segment.m:12
    segments = []
    # for i in arange(1,n).reshape(-1):
    for trace in f:
        # print(type(trace), trace)
        # clear('a','b','ya','yb','lin')
        a = np.nonzero(np.isfinite(trace.data))[0][0]
        b = np.nonzero(np.isfinite(trace.data))[0][-1]
        ya = trace.data[a]
        yb = trace.data[b]
        lin = np.linspace(start=ya, stop=yb, num=b-a+1)
        segment = trace.copy()
        # Ugly handling of different cases
        before = np.empty(a)
        before[:] = np.nan
        after = np.empty(len(trace.data)-b-1)
        after[:] = np.nan
        segment.data = np.concatenate((before, lin, after))
        segments.append(segment)
    if bare_trace:
        return segments[0]
    else:
        return segments


def _cum2grad(f_in, normalize=False, debug=False):
    """
    Transform cumulative into function that is shift at 0 for all maxima

    :param f_in: cumulative trace
    :param normalize: 'True' will divide the output by it's min value
    :returns: f_out
    """
    # cum2grad.m:7
    assert not len(f_in) == 0, 'f_in is empty!'

    tycalpha = np.zeros(len(f_in.data))
    tikxs, _ = _loca_ext(f_in, f_in.stats.starttime, f_in.stats.endtime, 'maxi')
    a = int(np.nonzero(np.isfinite(f_in.data))[0][0])
    b = int(np.nonzero(np.isfinite(f_in.data))[0][-1])
    # tikxs = tikxs.to_list()
    if len(tikxs) == 0:
        tikxs = np.array([a, b])
    else:
        tikxs = np.array([a] + tikxs.tolist() + [b])

    tikys = f_in.data[tikxs]
    # for j in np.arange(len(tikx) - 1, 1, -1).reshape(-1):
    for ilo, ihi, y in zip(tikxs[:-2], tikxs[1:], tikys[1:]):
        tycalpha[ilo:ihi] = y
        # tycalpha[np.arange(tikx(j), tikx(j + 1))] = tiky(j + 1)

    f_out = f_in.copy()
    f_out.data -= tycalpha
    f_out.data[f_out.data > 0] = 0
    if normalize:
        f_out.data /= abs(min(f_out.data))
    if debug:
        Stream([f_in, f_out]).plot()
    return f_out


def _loca_ext(trace, starttime, endtime, type_, debug=False):
    """
    Returns local extrema

    :param trace: waveform trace
    :param start_time: start of window to look at
    :param end_time: end of window to look at
    :param type_: 'maxi' or 'mini'
    :returns: [indice, value]
    """
    diff = trace.copy()
    diff.data = np.diff(np.sign(np.diff(trace.data)))
    diff = diff.slice(starttime, endtime)
    if debug:
        diff.plot()
    if type_ == 'maxi':
        loc = (diff.data < 0)
    else:
        loc = (diff.data > 0)

    i_extremes = np.nonzero(loc)[0] + 1
    return i_extremes, trace[i_extremes]


def _mean_trace(traces):
    """
    Calculate the mean trace

    :param traces: stream or list of traces
    :returns: trace object
    """
    if isinstance(traces, Trace):
        return traces.copy()
    if isinstance(traces[0], Trace):
        if len(traces) == 1:
            return traces[0].copy()
    else:
        assert isinstance(traces[0][0], Trace),\
            "traces not a Trace, list of Traces, or list of lists of Traces"
        traces = [t for x in traces for t in x]

    data_len = len(traces[0].data)
    for tr in traces[1:]:
        assert len(tr.data) == data_len, 'traces are not the same length'

    mean_tr = traces[0].copy()
    for tr in traces[1:]:
        mean_tr.data += tr.data
    mean_tr.data /= len(traces)
    return mean_tr


if __name__ == '__main__':
    pass
