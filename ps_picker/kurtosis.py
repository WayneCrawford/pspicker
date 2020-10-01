"""
Kurtosis-specific routines
"""
import warnings

import numpy as np
from scipy.signal import lfilter
from obspy.core.stream import Stream
from obspy.core.trace import Trace

from .trace_utils import same_inc, smooth_filter


class Kurtosis():
    """
    Kurtosis class

    englobes old trace2kurto, trace2FSkurto, fast_kurtosis f_cumul and
    smooth_filter.
    A bunch of static methods for now, until/less I find some common paramters
    """
    def __init__():
        pass

    @staticmethod
    def trace2kurto(trace, Fc, T, n_smooth):
        """
        Calculate cumulative kurtosis

        :param trace: waveform trace
        :param Fc: frequency bandwidth [low, high]
        :param T: sliding window length (s)
        :param n_smooth: smoothing sequence (samples)
        :returns: kurto_cum: cumulative kurtosis stream
                  ind_gmin: global minimum indices
                  gmin: global minimum value
        """
        sr = trace.stats.sampling_rate
        T_samples = np.floor(T * sr) + 1
        f = trace.copy()

        # Prep traces
        f.detrend('demean')
        f.filter(type='bandpass', freqmin=Fc[0], freqmax=Fc[1], corners=3)
        f.data[np.arange(0, 50)] = 0

        # Calculate the Kurtosis recursively
        kurtos = Kurtosis.fast_kurtosis(f, T_samples)
        f_kurto_cum = kurtos.copy()
        f_kurto_cum = smooth_filter(f_kurto_cum, n_smooth)
        kurto_cum, _ = Kurtosis.f_cumul(f_kurto_cum)
        line = Kurtosis.f_segment(kurto_cum)
        kurto_cum.data == line.data
        ind_gmin = kurto_cum.data.argmin()
        gmin = kurto_cum.data[ind_gmin]

        return kurto_cum, ind_gmin, gmin

    @staticmethod
    def trace2FWkurto(trace, FBs, window_lengths, n_smooth, starttime,
                      endtime):
        """
        Calculate cumulative kurtosis over bandwidths and windows

        FW stands for Frequency Window because it computes for
        several frequency bandwidth and window sizes

        :param trace: the raw trace (works one trace, not a stream)
        :param h: is the sampling frequency
        :param FBs: list of n [min, max] frequency bands
        :param window_lengths: list of m window lengths (seconds)
        :param starttime: is the first time of interest
        :param endtime: is the last time of interest
        :returns: mean kurtosis over all windows and frequency bandwidths
                  kurtosises for each window and frequency bandwidth
        """
        # trace2FWkurto.m:26
        sr = trace.stats.sampling_rate
        data_length = trace.stats.endtime - trace.stats.starttime
        if starttime < trace.stats.starttime:
            starttime < trace.stats.starttime
        if endtime > trace.stats.endtime:
            endtime = trace.stats.endtime

        # Filter traces in different frequency bands
        B = []
        for FB in FBs:
            f = trace.copy()
            f.filter(type='bandpass', freqmin=FB[0], freqmax=FB[1], corners=3)
            f.data = same_inc(f.data, starttime, endtime)
            B.append(f)

        K = []  # all kurtosises: 1st dim: window_lengths, 2nd dim: freq bands
        for win_len in window_lengths:
            if win_len > data_length:
                warnings.warn(
                    'Kurtosis window longer than data window '
                    '({:2g}s > {:4g}s), skipping!'.format(win_len,
                                                          data_length / sr))
            else:
                win_samps = np.floor(win_len * sr) + 1
                if win_samps > len(trace.data):
                    win_samps = len(trace.data)
                kurtos = [Kurtosis.fast_kurtosis(t, win_samps) for t in B]
                f_kurto_cums = smooth_filter(kurtos, n_smooth)
                kurto_cums = Kurtosis.f_cumul(f_kurto_cums)
                lines = Kurtosis.f_segment(kurto_cums)
                corr_kurto_cums = kurto_cums.copy()
                for c, l in zip(corr_kurto_cums, lines):
                    c -= l
                K.append(corr_kurto_cums)

        mean_K = K[0]
        mean_K.data = np.zeros(len(mean_K.data))
        for a in K:
            for b in a:
                mean_K.data += b.data
        mean_K.data /= (len(window_lengths) * len(FBs))
        return mean_K, K

    @staticmethod
    def fast_kurtosis(trace, win_samps):
        """
        Compute kurtosis really quickly using "filter" function

        could I just use obspy kurtosis?

        :param trace: one trace
        :param win_samps: number of samples in the sliding window
        :returns: Kurtosis trace
        """
        win_samps = int(round(win_samps))
        # fast_kurtosis.m:11
        if win_samps == 1:
            win_samps = 2

        # Set NaNs to 0 for computational stability
        f = trace.copy()
        f.detrend(type='demean')
        f.data[np.isnan(f.data)] = 0

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

    @staticmethod
    def follow_extrem(f, type_='mini', n_follow=2, smooth_vec=[1, 5, 10],
                      option=None, sense=None):
        """
        Find extremas using several smoothed versions of a function

        Steps:
          1) locate the 'n' first extremas of the smoothest function
          2) refine these points step by step through less and less smooth
             functions

        :param f: unsmoothed data trace
        :param type: 'mini' or 'maxi', depending on wheter you want to
            follow minima or maxima
        :param n_follow: number of extrema to follow
        :param smooth_vec: vector containing the different smoothings to apply
        :param option:
            'normalize': normalize the gradient
            otherwise:   don't normalize
        :param sense:
            'first': select first 'n_follow' extrema > 0.1, order right to left
            otherwise: select the 'n_follow' biggest extrema, order biggest
                       to smallest
        :returns:  smoothed Kurtosis traces (smoothest to roughest),
                   indices of extrema,
                   values of extrema
        """
        # follow_extrem.m:29
        # Parameters
        assert not len(f) == 0, 'f is empty!'
        assert isinstance(smooth_vec, list), "smooth_vec is not a list"

        # Holds a trace and the indices of its extrema for each smoothing
        smoothed = []
        smooth_vec.sort(reverse=True)  # Put the strongest smoothing first
        # for i in arange(1,length(smooth_vec)).reshape(-1):
        for smoothing in smooth_vec:
            v_smooth = smooth_filter(f, smoothing)
            v_smooth = Kurtosis.cum2grad(v_smooth, option)
            ext_indices, _ = Kurtosis._loca_ext(v_smooth, f.stats.starttime,
                                                f.stats.endtime, type_)
            smoothed.append({'indices': ext_indices, 'trace': v_smooth})

        # create a list of {index:  value:} dictionaries with all of the
        # detected extrema in the smoothest trace
        A = [{'index': i, 'value': smoothed[0]['trace'].data[i]}
             for i in smoothed[0]['indices']]
        # Choose up to 'n_follow' extrema in the smoothest trace
        if sense == 'first':
            # If extrema are tiny, take the biggest
            B = [x for x in A if x['value'] >= 0.1]  # eliminate tiny extrema?
            if len(B) == 0:
                iA = sorted(A, key=lambda k: np.abs(k['value']), reverse=True)
                try:
                    b = [x for x in iA[:n_follow]]
                except Exception:
                    b = [x for x in iA]
            else:
                iA = sorted(B, key=lambda k: k['index'])
                try:
                    b = [x for x in iA[:n_follow:-1]]
                except Exception:
                    b = [x for x in iA[::-1]]
        else:
            iA = sorted(A, key=lambda k: np.abs(k['value']), reverse=True)
            try:
                b = [x for x in iA[:n_follow]]
            except Exception:
                b = [x for x in iA]
        if len(b) == 0:
            return [x['trace'] for x in smoothed], b, []

        # Sharpen the indices/values using the smaller smoothing values
        # for k in arange(2, n).reshape(-1):
        for M, smoothing in zip(smoothed[1:], smooth_vec[1:]):
            if len(M['indices']) == 0:
                warnings.warn(
                    f'No extrema found for {smoothing:d}-sample '
                    'smoothing of Kurtosis')
            else:
                # for j in arange(1, len(b)).reshape(-1):
                c = []
                for extrem in b:
                    differ = abs(extrem['index'] - M['indices'])
                    if np.any(differ <= 40):
                        i = np.argmin(differ)
                        c.extend({'index': M['indices'][i],
                                  'value': M['trace'].data[i]})
                    else:
                        c.extend(extrem)
                b = c.copy()

        return ([x['trace'] for x in smoothed],
                [x['index'] for x in b],
                [x['value'] for x in b])

    @staticmethod
    def f_cumul(f):
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
            trace = t.copy()    # Backup copy for info
            t.data[np.isnan(t.data)] = 0
            t.differentiate(method='gradient')
            t.data[t.data < 0] = 0
            t.data = np.cumsum(t.data)
            t.data[np.isnan(trace)] = np.nan

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

    @staticmethod
    def f_segment(f):
        """
        Return a line segment between the first and last values of function.

        Goes from first non-nan value to last non-nan value
        Input and output have the same size.
        :param f: stream or list of data traces
        :returns: list of line segment traces
        """
        bare_trace = False
        if isinstance(f, Trace):
            f = [f]
            bare_trace = True
        # f_segment.m:12
        segments = []
        # for i in arange(1,n).reshape(-1):
        for trace in f:
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

    @staticmethod
    def cum2grad(f_in, option='', debug=False):
        """
        Transform cumulative into function that is shift at 0 for all maxima

        :param f_in: cumulative trace
        :param option: 'normalize' will divide the output by it's min value
        :returns: f_out
        """
        # cum2grad.m:7
        assert not len(f_in) == 0, 'f_in is empty!'

        tycalpha = np.zeros(len(f_in.data))
        tikxs, _ = Kurtosis._loca_ext(f_in, f_in.stats.starttime,
                                      f_in.stats.endtime, 'maxi')
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
        if option == 'normalize':
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
        diff = same_inc(diff, starttime, endtime)
        if debug:
            diff.plot()
        if type_ == 'maxi':
            loc = (diff.data < 0)
        else:
            loc = (diff.data > 0)

        i_extremes = np.nonzero(loc)[0] + 1
        return i_extremes, trace[i_extremes]


if __name__ == '__main__':
    pass
