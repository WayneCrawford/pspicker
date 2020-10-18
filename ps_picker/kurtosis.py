"""
Kurtosis-specific routines
"""
import warnings

import numpy as np
from scipy.signal import lfilter
from obspy.core.stream import Stream
from obspy.core.trace import Trace

from .trace_utils import smooth_filter
from .pick_candidate import PickCandidate


class Kurtosis_pick():
    """
    Holds one Kurtosis-based pick and the values that chose it
    """
    def __init__(self, time, k_jump, kurtosis):
        """
        :param time: kurtosis pick time
        :param k_jump: jump in kurtosis associated with the pick
        :param kurtosis: kurtosis value at pick
        """
        self.time = time
        self.k_jump = k_jump
        self.kurtosis = kurtosis


class Kurtosis():
    """
    Class for picking seismograms using the Kurtosis

    Kurtosis picking is a n-step process whose goal is to identify jumps
    in the Kurtosis:
    1) calculate the Kurtosis
    2) calculate the cumulative sum of the kurtosis (kurto_cum)
    3) remove a linear fit from the beginning to the end of the kurto_cum
       (corr_kurto_cum)
    4) Identify peaks (highs) and subtract this from the preceding data.  Set
       all positive values to zero
    5) identify where the slope of the resulting function changes from positive
       to negative
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

        # Mean trace over freq_bands, wind_lengths & smoothing
        self.mean_kurtosis = None
        # Mean trace, cumulated and detrended
        self.mean_cumulative_kurtosis = None
        # Gradients of mean_cumulative_kurtosis, for each extrem_smoothing
        self.kurto_gradients = None

    def __str__(self):
        s = "Kurtosis:\n"
        s += f"   freq bands = {self.freq_bands}\n"
        s += f"   window_lengths = {self.wind_lengths}\n"
        s += f"   n_smooth = {self.n_smooth}\n"
        s += f"   extrem_smooths = {self.extrem_smooths}\n"
        s += f"   mean_kurtosis = {self.mean_kurtosis}\n"
        s += f"   mean_cumulative kurtosis = {self.mean_cumulative_kurtosis}\n"
        s += f"   len{self.kurto_gradients} kurto_gradients"
        return s

    def pick_trace(self, trace, n_candidates, starttime=None, endtime=None,
                   extrem_type='mini', extrem_smooths=[1, 5, 10],
                   extrem_normalize=False, extrem_which='max'):
        """
        Pick a trace using the Kurtosis

        :param trace: waveform trace
        :param n_candidates: maximum number of candidates to return
        :param starttime: first time of interest
        :param endtime: last time of interest
        :param extrem_type: 'mini' or 'maxi', depending on wheter you want to
            follow minima or maxima
        :param extrem_smooths: smoothing to apply when following extrema
        :param extrem_normalize: normalize the cumulative Kurtosis gradient
        :param extrem_which:
            'first': select first 'n_follow' extrema > 0.1, order right to left
            'max': select the 'n_follow' biggest extrema, order biggest
                       to smallest
        :returns:  list of PickCandidates
        """
        self.extrem_smooths = extrem_smooths
        self.calc_kurtocum(trace, starttime, endtime)
        candidates = self.follow_extrem(type_=extrem_type,
                                        n_follow=n_candidates,
                                        normalize=extrem_normalize,
                                        sense=extrem_which)
        return candidates

    def calc_kurtocum(self, trace, starttime=None, endtime=None, debug=False):
        """
        Calculate cumulative kurtosis over bandwidths and windows

        Puts mean kurtosis over all windows and freq_bands in
        self.mean_cumulative_kurtosis
        :param trace: the raw trace (works one trace, not a stream)
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
        Return extrema of the cumulative kurtosis

        Steps:
          1) Smooth the cumulative kurtosis (self.mean_cumulative_kurtosis)
             using the windows specified in self.extrem_smooths
          2) locate the first 'n' first extremas of the smoothest function
          3) refine these extrema step by step through less and less smooth
             functions

        Uses self.mean_cumulative kurtosis
        :param type: 'mini' or 'maxi':  follow minima or maxima
        :param n_follow: number of extrema to follow
        :param normalize: normalize the gradient?
        :param sense:
            'first': select first 'n_follow' extrema > 0.1, ordered right
                     to left
            'max': select from biggest to smallest
        :returns:  list of PickCandidates
        """
        # follow_extrem.m:29
        # Parameters
        assert not len(self.mean_cumulative_kurtosis) == 0,\
            'no mean cumulative kurtosis!'
        st = self.mean_cumulative_kurtosis.stats.starttime
        sr = self.mean_cumulative_kurtosis.stats.sampling_rate

        all_extrema, self.kurto_gradients, _ =\
            _get_extrema(self.extrem_smooths, self.mean_cumulative_kurtosis,
                         type_, normalize)

        selected_extrema = _select_extrema(all_extrema[0],
                                           n_follow, sense)

        # Sharpen the indices/values using the smaller smoothing values
        sharp_extrema = []
        for e in selected_extrema:
            extrema = e
            for finer_extrema in all_extrema[1:]:
                c_extrem = _find_close_extrema(
                    extrema, finer_extrema, 40)
                if c_extrem is not None:
                    extrema = c_extrem
            sharp_extrema.append(extrema)

        return [PickCandidate(st + x['index']/sr,
                              'kurtosis',
                              x['value'],
                              sampling_rate=sr)
                for x in sharp_extrema]


def _find_close_extrema(best_extrem, finer_extrema, max_diff=40):
    if len(finer_extrema) == 0:
        warnings.warn('No extrema found for {:d}-sample smoothing'.format(
                      finer_extrema['smoothing']))
        return None
    i_diff = abs(best_extrem['index'] - [x['index'] for x in finer_extrema])
    if np.any(i_diff <= max_diff):
        return finer_extrema[np.argmin(i_diff)]
    return None


def _get_extrema(smoothing_list, cumul_k, sense, normalize, debug=False):
    """
    Returns a list of extrema, from smoothest to roughest

    :param smoothing list: list of smoothings in number of samples
    :param cumul_k: cumulative kurtosis Trace
    :param sens: 'maxi' or 'mini', passed on to ext_indices
    :param normalize: passed on to _cum2grad()
    :returns: list of {indices: value, trace: value} from smoothest to
        roughest.  value is the approximate height of the kurtosis jump
    """
    extrema, gradients = [], []
    # Put the strongest smoothing first
    sorted_smooth = sorted(smoothing_list, reverse=True)
    for smoothing in sorted_smooth:
        v_smooth = smooth_filter(cumul_k, smoothing)
        v_smooth = _cum2grad(v_smooth, normalize)
        ext_indices, _ = _loca_ext(v_smooth, None, None, sense)
        extrema.append([{'index': i, 'value': -v_smooth.data[i]}
                        for i in ext_indices])
        gradients.append(v_smooth)
    # create a list of {index:  value:} dictionaries with all of the
    # detected extrema in the smoothest trace
    # extrema = [{'index': i, 'value': extrema[0]['trace'].data[i]}
    #            for i in extrema[0]['indices']]
    return extrema, gradients, sorted_smooth


def _select_extrema(extrema, N, sense='max', threshold=0.1):
    """
    Return N extrema, ordered by size or time
    :param extrema: list of {'value': val, 'index': i}
    :param sense: 'max' return from largest to smallest
                  'first': return from first to last
    :param threshold: minimum value to accept for sense=='first'
    """
    if sense == 'first':
        big_extrema = [x for x in extrema if x['value'] >= threshold]
        if len(big_extrema) > 1:
            ext = sorted(big_extrema, key=lambda k: k['index'])
        else:
            ext = sorted(extrema, key=lambda k: k['index'])
        try:
            selected = [x for x in ext[:N:-1]]
        except Exception:
            selected = [x for x in ext[::-1]]
    elif sense == 'max':
        ext = sorted(extrema, key=lambda k: np.abs(k['value']), reverse=True)
        try:
            selected = [x for x in ext[:N]]
        except Exception:
            selected = [x for x in ext]
    else:
        raise NameError("sense not 'max' or 'first'")
    return selected


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
    # an NaN in the original data.
    # I think this is outdated, should just trim
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
    Transform cumulative into function that is shifted below 0 for all maxima

    Finds peaks, subtracts the next peak value from all data, gets rid of
    all values > 0.  Give negative values that "jump" up to zero at the next
    peak

    :param f_in: cumulative trace
    :param normalize: 'True' will divide the output by it's min value
    :returns: f_out
    """
    # cum2grad.m:7
    assert not len(f_in) == 0, 'f_in is empty!'

    tycalpha = np.zeros(len(f_in.data))
    tikxs, _ = _loca_ext(f_in, None, None, 'maxi')
    a = int(np.nonzero(np.isfinite(f_in.data))[0][0])
    b = int(np.nonzero(np.isfinite(f_in.data))[0][-1])
    # tikxs = tikxs.to_list()
    if len(tikxs) == 0:
        tikxs = np.array([a, b])
    else:
        tikxs = np.array([a] + tikxs.tolist() + [b])

    tikys = f_in.data[tikxs]
    # function equal to the next peak
    for ilo, ihi, y in zip(tikxs[:-2], tikxs[1:], tikys[1:]):
        tycalpha[ilo:ihi] = y
        # tycalpha[np.arange(tikx(j), tikx(j + 1))] = tiky(j + 1)

    f_out = f_in.copy()
    f_out.data -= tycalpha          # input minus the next peak value
    f_out.data[f_out.data > 0] = 0  # Get rid of everything above zero
    if normalize:
        f_out.data /= abs(min(f_out.data))
    if debug:
        Stream([f_in, f_out]).plot()
    return f_out


def _loca_ext(trace, starttime, endtime, type_, debug=False):
    """
    Returns local extrema

    Actually just returns where the trace slope changes from positive to
    negative (type_=='maxi'), or vice versa (type='mini')
    :param trace: waveform trace
    :param start_time: start of window to look at
    :param end_time: end of window to look at
    :param type_: 'maxi' or 'mini'
    :returns: indice, trace_value, diff_value
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
    return i_extremes, trace.data[i_extremes]


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
