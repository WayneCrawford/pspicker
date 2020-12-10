import numpy as np
from scipy import signal, linalg
from obspy.core import UTCDateTime
from obspy.core.stream import Stream

from .logger import log
# from .timer import Timer


class Polarity:
    """
    polarity and dip-rect calculation
    """
    def __init__(self, stream, params, zcomponents='Z3', ncomponents='N1Y',
                 ecomponents='E2X', verbose=True):
        """
        :param stream: stream containing X, Y and Z traces
        :param params: PolarityParameters object
        :param zcomponents: possible characters for vertical components
        :param ncomponents: possible characters for "north" components
        :param ecomponents: possible characters for "east" components
        :param
        """
        stream, self.nsamps = self._clean_stream(stream)
        self._check_components(zcomponents, ncomponents, ecomponents)
        self.params = params
        self.verbose = verbose
        self.sr = stream[0].stats.sampling_rate
        self.tracez = None
        self.tracen = None
        self.tracee = None
        for tr in stream:
            if tr.stats.channel[-1] in zcomponents:
                assert self.tracez is None, 'already have a z trace!'
                self.tracez = tr.copy()
            elif tr.stats.channel[-1] in ncomponents:
                assert self.tracen is None, 'already have a n trace!'
                self.tracen = tr.copy()
            elif tr.stats.channel[-1] in ecomponents:
                assert self.tracee is None, 'already have a e trace!'
                self.tracee = tr.copy()
            else:
                msg = f"Unexpected channel code {tr.stats.channel}"
                raise ValueError(msg)

    @staticmethod
    def _clean_stream(st):
        """
        Validate stream and synchronize traces
        """
        assert len(st) == 3, 'stream does not have three channels'
        assert st[0].stats.sampling_rate == st[1].stats.sampling_rate,\
            'traces have different sampling rates'
        assert st[1].stats.sampling_rate == st[2].stats.sampling_rate,\
            'traces have different sampling rates'

        stime = UTCDateTime(np.max([x.stats.starttime.timestamp for x in st]))
        etime = UTCDateTime(np.min([x.stats.endtime.timestamp for x in st]))
        st.trim(stime, etime)
        n_samps = np.min([len(x.data) for x in st])
        for tr in st:
            if not len(tr.data) == n_samps:
                tr.data = tr.data[:n_samps]
        return st, n_samps

    @staticmethod
    def _check_components(zc, nc, ec):
        """
        Verify that there are no repeated characters in the component maps

        :param zc: characters allowed for z component
        :param nc: characters allowed for n component
        :param ec: characters allowed for e component
        """
        error_str = "overlap between {} and {} component maps: {}, {}"
        for c in zc:
            assert c not in nc, error_str.format('z', 'n', zc, nc)
            assert c not in ec, error_str.format('z', 'e', zc, ec)
        for c in ec:
            assert c not in nc, error_str.format('e', 'n', ec, nc)

    def calc_dip_rect(self, times, min_dip_threshold=30, plot=False):
        """
        Return dip-rect values after specified times

        Dip-rect is the rectilinearity multiplied by -1 for near-horizontal
            particle motions, by 1 for near-vertical particle motions.

        :param times: list of UTCDateTimes
        :param min_dip_threshold: do not calculate dip-rectilinearity if there
            is no calculated dip with at least this absolute angle (degrees)
        :returns: DR: Dip-rectilinearity trace
        """
        if len(times) == 0:
            return None
        # Pick_Function.m:543
        # with Timer(text="    polarity.calc_dip_rect(): {:0.4f}s"):
        rectP, aziP, dipP = self.polar_analysis(times)
        if np.max(np.abs(dipP.data)) < min_dip_threshold:
            if self.verbose:
                log("max dip ({:.1f}) < threshold ({:g}), no DR calculated".
                    format(np.max(np.abs(dipP.data)), min_dip_threshold),
                    "verbose")
            return None

        dipp = np.sin(np.abs(np.radians(dipP.data)))  # 1 for vert, 0 for hor
        smooth_samples = int(np.round(self.params.DR_smooth_length * self.sr))
        if smooth_samples < 2:
            smooth_samples = 2
        a = np.ones(smooth_samples) / smooth_samples
        a2 = np.ones(2*smooth_samples) / (2 * smooth_samples)
        smooth_dipp = signal.lfilter(a, 1, dipp)
        smooth_rectP = signal.lfilter(a, 1, rectP.data)
        Drb = np.sign(1.3 * smooth_dipp - smooth_rectP)  # pos for P, neg for S
        DR = rectP.copy()
        DR.data *= Drb
        DR.data = signal.lfilter(a2, 1, DR.data)
        # DR.data = signal.lfilter(a2, 1, np.multiply(rectP.data, Drb))
        DR.data[rectP.data == 0] = 0
        if plot:
            DR.stats.channel = 'DR'
            dipP.stats.channel = 'DIP'
            Stream([DR, dipP, self.tracez, self.tracen, self.tracee]).plot(
                equal_scale=False)
        return DR

    def polar_analysis(self, times):
        """
        Compute the polarity parameters around specified times

        Only computes around user specified times since this analysis is time
        consuming.  Other values are set to zero
        :param times: list of times around which to compute polarity
        :returns: rect, azi (degrees), dip (degrees)
        """
        # fast_polar_analysis.m:17
        # with Timer(text="polarity.polar_analysis() calc_indices: {:0.4f}s"):
        ind_vec, n_half_analyze = self._calc_indices(times)

        # with Timer(text="    polarity.polar_analysis() rest: {:0.4f}s"):
        rectP = self.tracez.copy()
        rectP.data[:] = 0.
        aziP = rectP.copy()
        dipP = rectP.copy()
        for k in ind_vec:
            eP = self.tracee.data[k - n_half_analyze: k + n_half_analyze]
            nP = self.tracen.data[k - n_half_analyze: k + n_half_analyze]
            zP = self.tracez.data[k - n_half_analyze: k + n_half_analyze]
            MP = np.cov(np.array([eP, nP, zP]))
            D, V = linalg.eig(MP)
            i_sort = np.argsort(D)  # sort from min (0) to max (2)
            D = np.abs(D[i_sort])
            V = V[:, i_sort]
            # 2 = major axis, 0 = minor axis, 1 = intermediate
            rectP.data[k] = 1 - ((D[0] + D[1]) / (2 * D[2]))
            V_major = V[:, 2]
            aziP.data[k] = np.degrees(np.arctan(V_major[1] / V_major[0]))
            dipP.data[k] = np.degrees(np.arctan(V_major[2]
                                                / np.sqrt(V_major[1]**2
                                                          + V_major[0]**2)))
        return rectP, aziP, dipP

    def _calc_indices(self, pick_times):
        """
        Create list of indices to investigate

        Compute windows after picks minus duplicate samples or samples too
        close to the data edges
        """
        n_compute = round(self.sr * self.params.calculate_window)
        ind_vec = np.array([], dtype='int32')
        # ztimes = self.tracez.times('utcdatetime')
        ztimes = self.tracez.times('timestamp')
        for t in pick_times:
            index = ztimes.searchsorted(t.timestamp)
            iwind = np.arange(index, index + n_compute, dtype='int32')
            ind_vec = np.concatenate([ind_vec, iwind])

        n_half_analyze = int(round(self.params.analyze_window * self.sr / 2))
        # Get rid of duplicate or out-of-bounds values
        ind_vec = ind_vec[ind_vec - n_half_analyze >= 0]
        ind_vec = ind_vec[ind_vec + n_half_analyze < self.nsamps]
        ind_vec = np.sort(ind_vec)

        # ind_vec[np.diff(ind_vec, prepend=-1) == 0] = []
        ind_vec = ind_vec[np.diff(ind_vec, prepend=-1) != 0]
        # print(ind_vec, n_half_analyze)
        return ind_vec, n_half_analyze

    if __name__ == '__main__':
        pass
