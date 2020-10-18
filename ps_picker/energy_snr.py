import copy

import numpy as np
from obspy.signal.util import smooth as obspy_smooth
from obspy.core.stream import Stream

from .trace_utils import smooth_filter


class EnergySNR():
    """
    Energy and signal-to-noise ratio of a waveform Stream
    """
    def __init__(self, stream, params, plot=False):
        """
        calculate energy and signal-to-noise ratio of a seismic stream

        :param stream: group of related traces to calculate over
        :param params: SNRParameters object
        """
        assert isinstance(stream, Stream), "stream is not an obspy Stream"

        self.params = params
        self.snr_threshold = None
        self.nrg = self._calc_energy(stream)
        self.snr = self._calc_snr()
        self.nrg.stats.channel = 'NRG'
        self.snr.stats.channel = 'SNR'

        if plot:
            snr_dB = self.snr.copy()
            snr_dB.data = 20*np.log10(snr_dB.data)
            snr_dB.stats.channel = 'SDB'
            Stream([snr_dB, self.snr, self.nrg]).plot(equal_scale=False)

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)

    @staticmethod
    def _calc_energy(stream):
        temp = stream.copy()
        for t in temp:
            t.data = np.power(t.data, 2)
        energy = temp.stack()[0]
        energy.data = np.power(energy.data, 0.5)
        return energy

    def _calc_snr(self):
        """
        Return estimated signal to noise ratio

        :param energy: energy trace
        :returns: (M(t...t+window_after)/M(t-window_before,t))
        """
        sr = self.nrg.stats.sampling_rate
        half_noise_wind_samps = round(sr * self.params.noise_window / 2)
        half_signal_wind_samps = round(sr * self.params.signal_window / 2)
        self.params.noise_window = 2*half_noise_wind_samps/sr
        self.params.signal_window = 2*half_signal_wind_samps/sr
        w_noise = self.nrg.copy()
        w_signal = self.nrg.copy()
        # smooth using central moving average
        w_noise.data = obspy_smooth(self.nrg.data, half_noise_wind_samps)
        w_signal.data = obspy_smooth(self.nrg.data, half_signal_wind_samps)
        # Shift times so that noise is BEFORE reference time and signal AFTER
        w_noise.stats.starttime += self.params.noise_window/2
        w_signal.stats.starttime -= self.params.signal_window/2
        global_start = w_noise.stats.starttime
        global_end = w_signal.stats.endtime
        w_noise.trim(global_start, global_end)
        w_signal.trim(global_start, global_end)
        w_signal.data /= w_noise.data
        return w_signal

    def slice(self, starttime=None, endtime=None, nearest_sample=True):
        """
        Returns sliced view of the nrg and snr traces
        """
        new = self.copy()
        new.nrg = self.nrg.slice(starttime=starttime, endtime=endtime,
                                 nearest_sample=nearest_sample)
        new.snr = self.snr.slice(starttime=starttime, endtime=endtime,
                                 nearest_sample=nearest_sample)
        return new

    def is_trustworthy(self, n_smooth=100, debug=False):
        """
        Check if the signal can be trusted or not.

        Consider the signal trustworthy if the SNR crosses (from below to
        above) the specified threshold at least once and no more than
        SNR_threshold_crossings times
        :param n_smooth: length of moving average filter to apply before
            analysis
        """
        # Pick_Function.m:380
        snr_smooth = smooth_filter(self.snr, n_smooth)
        self._set_snr_threshold(snr_smooth)
        sign_change = np.diff(np.sign(snr_smooth.data - self.snr_threshold))
        crossings = len(sign_change[sign_change == 2])
        if debug:
            snr_smooth.plot()
        return (crossings > 0
                and crossings <= self.params.max_threshold_crossings)

    def _set_snr_threshold(self, snr_smooth):
        """
        Calculate the signal-to-noise threshold value

        :param snr_smooth: trace of smoothed signal-to-noise ratio
        """
        # Pick_Function.m:382
        tp = self.params.threshold_parameter
        min_threshold = min(self.params.quality_thresholds)
        if (tp > 0 and tp <= 1):
            threshold = 1 + tp * (np.nanmax(snr_smooth.data) - 1)
        else:
            assert tp < 0, f'Illegal SNR threshold_parameter value: {tp:g}'
            threshold = -tp
        self.snr_threshold = max([threshold, min_threshold])


if __name__ == '__main__':
    pass
