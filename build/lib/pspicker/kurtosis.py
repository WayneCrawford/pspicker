"""
Kurtosis-specific routines
"""
import warnings

import numpy as np
from scipy.signal import lfilter
from obspy.core.stream import Stream
from obspy.core.trace import Trace
# from scipy.stats import kurtosis as scipy_kurtosis
# from obspy.realtime.signal import kurtosis as obspy_kurtosis

from .utils import smooth_filter
from .pick_candidate import PickCandidate
from .logger import log


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
    def __init__(self, params, plot=False):
        """
        :param params: KurtosisParameters object
        :param plot: make a plot of Kurtosis parameters
        """
        self.params = params
        self.plot = plot

        # Mean trace over freq_bands, wind_lengths & smoothing
        self.mean_kurtosis = None
        # Mean trace, cumulated and detrended
        self.mean_cumulative_kurtosis = None
        # Gradients of mean_cumulative_kurtosis, for each extrem_smoothing
        self.kurto_gradients = None

    def __str__(self):
        s = "Kurtosis:\n"
        s += f"   params = {self.params}\n"
        s += f"   mean_kurtosis = {self.mean_kurtosis}\n"
        s += f"   mean_cumulative kurtosis = {self.mean_cumulative_kurtosis}\n"
        s += f"   len{self.kurto_gradients} kurto_gradients"
        return s

    def pick_trace(self, trace, n_candidates, starttime=None, endtime=None,
                   extrem_type='mini', extrem_normalize=False,
                   extrem_which='max'):
        """
        Pick a trace using the Kurtosis

        :param trace: waveform trace
        :param n_candidates: maximum number of candidates to return
        :param starttime: first time of interest
        :param endtime: last time of interest
        :param extrem_type: 'mini' or 'maxi', depending on wheter you want to
            follow minima or maxima
        :param extrem_normalize: normalize the cumulative Kurtosis gradient
        :param extrem_which:
            'first': select first 'max_cadidates' extrema > 0.1, order right to
                     left
            'max': select the 'max_candidates' biggest extrema, order biggest
                       to smallest
        :returns:  list of PickCandidates
        """
        self.calc_kurtocum(trace, starttime, endtime)
        candidates = self.follow_extrem(ext_type=extrem_type,
                                        max_candidates=n_candidates,
                                        normalize=extrem_normalize,
                                        order=extrem_which)
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

        # sr = trace.stats.sampling_rate
        # data_length = trace.stats.endtime - trace.stats.starttime
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
        log('Pre-filtering data for kurtosis in {} bands'.format(
            self.params.frequency_bands), 'debug')
        for FB in self.params.frequency_bands:
            f = trace.copy()
            f.detrend('demean')
            # RAMP UP SIGNAL BEFORE FILTERING
            # ramp_time = 2 / FB[0]
            # ramp_pts = int(max(ramp_time * f.stats.sampling_rate, 50))
            # ramp = np.arange(0, ramp_pts) / ramp_pts
            # f.data[:ramp_pts] *= ramp
            f.filter(type='bandpass', freqmin=FB[0], freqmax=FB[1], corners=3)
            f = f.slice(starttime, endtime)
            # RAMP UP SIGNAL AFTER FILTERING
            # f.data[:ramp_pts] *= ramp
            B.append(f)

        # 2-level lists: : 1st dim: window_lengths, 2nd dim: freq bands
        K, C = [], []  # all kurtoses and all cumulative, detrended kurtoses
        log('Calculating kurtosis of filtered data using {}-s windows'
            .format(self.params.window_lengths), 'debug')
        for win_len in self.params.window_lengths:
            corr_cums, kurtos = self.calc_cum_kurtoses(B, win_len)
            C.extend(corr_cums)
            K.extend(kurtos)
        self.mean_kurtosis = _mean_trace(K)
        self.mean_cumulative_kurtosis = _mean_trace(C)

    def calc_cum_kurtoses(self, B, win_len):
        """
        Calculate kurtoses and cumulative kurtoses for given window length

        :param B: list of traces, filtered in different bands
        :win_len: window length in seconds
        """
        corr_cums, kurtos = [], []
        for tr, fb in zip(B, self.params.frequency_bands):
            sr = tr.stats.sampling_rate
            dl = tr.stats.endtime - tr.stats.starttime
            if win_len > dl:
                warnings.warn('Kurtosis window > data window ('
                              f'{win_len:.3g}s > {dl:.3g}s), skipping!')
                continue
            win_samps = min(int(np.floor(win_len * sr)) + 1, len(tr.data))
            k = _fast_kurtosis(tr, win_samps)
            filt = smooth_filter(k, self.params.n_smooth)
            corr_cum, _ = _f_cumul(filt.copy())
            corr_cum.detrend('simple')
            if self.plot:
                self._plot_kurtosis(win_len, fb, tr, filt, corr_cum)
            kurtos.append(k)
            corr_cums.append(corr_cum)
        return corr_cums, kurtos

    @staticmethod
    def _plot_kurtosis(wl, fb, trace, kurtosis, corr_cum):
        """
        :param wl: window length (s)
        :param fb: frequency band
        :param trace: data trace
        :param kurtosis: filtered Kurtosis
        :param corr_cum: corrected cumulative Kurtosis
        """
        log(f'Plot kurtosis for win_len={wl}, f_band={fb}', 'debug')
        cor = corr_cum.copy()
        cor.stats.channel = 'COR'
        kurtosis.stats.channel = 'KUR'
        Stream([trace, kurtosis, cor]).plot(
            size=(600, 600), equal_scale=False)

    def follow_extrem(self, ext_type='mini', max_candidates=2,
                      normalize=False, order=None, debug=False):
        """
        Return extrema of the cumulative kurtosis

        Steps:
          1) Smooth the cumulative kurtosis (self.mean_cumulative_kurtosis)
             using the windows specified in self.params.extrema_smoothings
          2) locate the first 'n' first extremas of the smoothest function
          3) refine these extrema step by step through less and less smooth
             functions

        Uses self.mean_cumulative kurtosis
        :param ext_type: 'mini' or 'maxi':  follow minima or maxima
        :param max_candidates: max number of candidates to return
        :param normalize: normalize the gradient?
        :param order:
            'first': select first 'max_candidates' extrema > 0.1, ordered right
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

        all_extrema = self._get_extrema(ext_type, normalize)

        selected_extrema = _select_extrema(all_extrema[0],
                                           max_candidates,
                                           order)

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

    def _get_extrema(self, ext_type, normalize, debug=False):
        """
        Reverse sorts extrema_smoothings and returns corresponding extrema

        :param ext_type: 'maxi' or 'mini', passed on to ext_indices
        :param normalize: passed on to _cum2grad()
        :returns: list of lists of {indices: value, trace: value} from
            smoothest to roughest.  value is approx height of the kurtosis jump
        """
        extrema, gradients = [], []
        # Put the strongest smoothing first
        self.params.extrema_smoothings.sort(reverse=True)
        # sorted_smooth = sorted(self.params.extrema_smoothings, reverse=True)
        # self.params.extrema_smoothings = sorted_smooth = sorted(
        for smoothing in self.params.extrema_smoothings:
            v_smooth = smooth_filter(self.mean_cumulative_kurtosis, smoothing)
            v_smooth = _cum2grad(v_smooth, normalize)
            ext_indices = _loca_ext(v_smooth, ext_type)
            extrema.append([{'index': i, 'value': -v_smooth.data[i]}
                            for i in ext_indices])
            gradients.append(v_smooth)
        self.kurto_gradients = gradients
        return extrema


def _find_close_extrema(best_extrem, finer_extrema, max_diff=40):
    if len(finer_extrema) == 0:
        warnings.warn('No extrema found for {:d}-sample smoothing'.format(
                      finer_extrema['smoothing']))
        return None
    i_diff = abs(best_extrem['index'] - [x['index'] for x in finer_extrema])
    if np.any(i_diff <= max_diff):
        return finer_extrema[np.argmin(i_diff)]
    return None


def _select_extrema(extrema, N, order='max', threshold=0.1):
    """
    Return N extrema, ordered by size or time
    :param extrema: list of {'value': val, 'index': i}
    :param order: 'max' return from largest to smallest
                  'first': return from first to last
    :param threshold: minimum value to accept for ex=='first'
    """
    if order == 'first':
        big_extrema = [x for x in extrema if x['value'] >= threshold]
        if len(big_extrema) > 1:
            ext = sorted(big_extrema, key=lambda k: k['index'])
        else:
            ext = sorted(extrema, key=lambda k: k['index'])
        try:
            selected = [x for x in ext[:N:-1]]
        except Exception:
            selected = [x for x in ext[::-1]]
    elif order == 'max':
        ext = sorted(extrema, key=lambda k: np.abs(k['value']), reverse=True)
        try:
            selected = [x for x in ext[:N]]
        except Exception:
            selected = [x for x in ext]
    else:
        raise NameError("order not 'max' or 'first'")
    return selected


def _fast_kurtosis(trace, win_samps):
    """
    Compute kurtosis quickly using "filter" function

    could I just use obspy kurtosis?

    :param trace: one trace
    :param win_samps: number of samples in the sliding window
    :returns: Kurtosis trace
    """
    assert isinstance(trace, Trace), "trace is not an obspy Trace"
    out = trace.copy()
    win_samps = int(round(win_samps))
    # log(f'{win_samps=}', 'debug')
    if win_samps == 1:
        win_samps = 2
    a = np.divide(np.ones(win_samps), float(win_samps))
    b = 1.

    f = trace.copy()
    f.detrend(type='demean')
    # Make buffer using first second of data
    one_sec = int(trace.stats.sampling_rate)
    buffer = np.tile(f.data[:one_sec], int(np.ceil(win_samps/one_sec)))
    data = np.concatenate((buffer[:win_samps], f.data))
    # data = np.concatenate((np.ones(win_samps)*f.data[0], f.data))

    # Compute kurtosis
    m_2 = lfilter(a, b, data**2)
    m_4 = lfilter(a, b, data**4)
    out.data = np.divide(m_4, (m_2 ** 2))
    # Cut off buffer
    out.data = out.data[win_samps:]
    # Protect against edge effect
    # out.data[:win_samps] = out.data[win_samps]

    # Set any kurtosis value to nan for any indices within win_samples of
    # an NaN in the original data.
    # I think this is outdated, should just trim
    # for i in np.nonzero(trace.data == np.nan)[0]:
    #     if i:
    #         out.data[i: i + win_samps] = np.nan
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
    :param f: trace, stream or list of traces containing the data
    :returns: trace or list of cumulative output traces, first value = 0
              list of cumulative output traces
    """
    # f_cumul.m:14
    bare_trace = False
    if isinstance(f, Trace):
        f = [f]
        bare_trace = True
    g = f.copy()
    for t in g:
        # tdata = t.data.copy()    # Backup copy for info
        t.data[np.isnan(t.data)] = 0
        t.data = np.gradient(t.data)
        # t.data = np.diff(t.data, prepend=0)
        # t.differentiate(method='gradient')
        t.data[t.data < 0] = 0
        t.data = np.cumsum(t.data)
        # t.data[np.isnan(tdata)] = np.nan

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
    :param f: trace, stream or list of data traces
    :returns: line segment trace or list of traces
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
    tikxs = _loca_ext(f_in, 'maxi')
    # a = int(np.nonzero(np.isfinite(f_in.data))[0][0])
    # b = int(np.nonzero(np.isfinite(f_in.data))[0][-1])
    a = 0
    b = len(f_in.data) - 1
    # print(f'{f_in.data.shape=}, {tycalpha.shape=}, {b=}')
    # print(f'{tikxs=}')
    # tikxs = tikxs.to_list()
    if len(tikxs) == 0:
        tikxs = np.array([a, b])
    else:
        tikxs = np.array([a] + tikxs.tolist() + [b])

    # function equal to the next peak
    for j in range(len(tikxs)-2, -1, -1):
        tycalpha[tikxs[j]:tikxs[j+1]+1] = f_in.data[tikxs[j+1]]

    f_out = f_in.copy()
    f_out.data -= tycalpha          # input minus the next peak value
    f_out.data[f_out.data > 0] = 0  # Get rid of everything above zero
    if normalize:
        f_out.data /= abs(min(f_out.data))
    if debug:
        Stream([f_in, f_out]).plot()
    return f_out


def _loca_ext(trace, ext_type, starttime=None, endtime=None, debug=False):
    """
    Returns local extrema

    Actually just returns where the trace slope changes from positive to
    negative (ext_type=='maxi'), or vice versa (type='mini')
    :param trace: waveform trace
    :param ext_type: 'maxi' or 'mini'
    :param start_time: start of window to look at
    :param end_time: end of window to look at
    :returns: indices
    """
    assert ext_type in ('mini', 'maxi')
    diff = trace.copy()
    diff.data = np.diff(np.sign(np.diff(trace.data)), prepend=0)
    if starttime is not None or endtime is not None:
        diff.trim(starttime, endtime)
    if debug:
        diff.plot()
    if ext_type == 'maxi':
        loc = (diff.data < 0)
    else:
        loc = (diff.data > 0)

    i_extremes = np.nonzero(loc)[0] + 1
    return i_extremes


def _mean_trace(traces):
    """
    Calculate the mean trace

    :param traces: stream or list of traces
    :returns: trace object
    """
    if isinstance(traces, Trace):
        return traces.copy()
    elif isinstance(traces[0], Trace):
        if len(traces) == 1:
            return traces[0].copy()
    else:
        ValueError("traces is not a Trace or a list of Traces")

    s = Stream(traces).copy()
    s = s.stack(npts_tol=1, time_tol=1./traces[0].stats.sampling_rate)
    if s[0].stats.starttime.timestamp == 0:
        ValueError('input traces did not start at same time')
    return s[0]
    # data_len = len(traces[0].data)
    # for tr in traces[1:]:
    #     assert len(tr.data) == data_len, 'traces are not the same length'
    #
    # log(f'Calculating mean of {len(traces):d} traces', 'debug')
    # mean_tr = traces[0].copy()
    # for tr in traces[1:]:
    #     mean_tr.data += tr.data
    # mean_tr.data /= len(traces)
    # return mean_tr


if __name__ == '__main__':
    pass
