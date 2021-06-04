import copy

import numpy as np
from obspy.signal.util import smooth as obspy_smooth
from obspy.core.stream import Stream
from obspy.core import UTCDateTime

from .utils import smooth_filter
from .logger import log


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
        self.station = self._get_station(stream)
        self.snr = self._calc_snr()
        self.nrg.stats.channel = 'NRG'
        self.snr.stats.channel = 'SNR'

        if plot:
            snr_dB = self.snr.copy()
            snr_dB.data = 20*np.log10(snr_dB.data)
            snr_dB.stats.channel = 'SDB'
            Stream([stream[0], snr_dB, self.snr, self.nrg]).plot(
                equal_scale=False)

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)

    @staticmethod
    def _calc_energy(stream):
        temp = stream.copy()
        for t in temp:
            t.data = np.power(t.data, 2)
        try:
            energy = temp.stack(npts_tol=1,
                                time_tol=1/temp[0].stats.sampling_rate)[0]
            if energy.stats.starttime.timestamp == 0:
                log('stacked traces are offset by more than one sample'
                    'setting starttime to stream[0].stats.starttime',
                    'error')
                energy.stats.starttime = stream[0].stats.starttime
        except ValueError as err:
            log(err, 'error')
            log('Slicing to same size', 'warning')
            log(temp, 'warning')
            last_start = UTCDateTime(max([t.stats.starttime.timestamp
                                          for t in temp]))
            first_end = UTCDateTime(min([t.stats.endtime.timestamp
                                         for t in temp]))
            temp.trim(last_start, first_end)
            log(temp, 'warning')
            energy = temp.stack(npts_tol=1)[0]
        energy.data = np.power(energy.data, 0.5)
        return energy

    @staticmethod
    def _get_station(stream):
        stations = list(set([t.stats.station for t in stream]))
        if len(stations) == 1:
            return stations[0]
        else:
            return "MULTIPLE STATIONS"

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
        noise_start = w_noise.stats.starttime
        signal_end = w_signal.stats.endtime
        if noise_start >= signal_end:
            log("SNR noise_start >= signal_end, no SNR calculated"
                "(trace, signal, noise lengths = {.2g}, {.2g}, {.2g})".format(
                self.nrg.stats.npts / self.nrg.stats.sampling_rate,
                self.params.signal_window,
                self.params.noise_window), 'error')
            return None
        w_noise.trim(noise_start, signal_end)
        w_signal.trim(noise_start, signal_end)
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
        :returns: is_trustworthy, info_text
        """
        # Pick_Function.m:380
        if self.snr.stats.npts == 0:
            log(f'station {self.station} self.snr has zero length', 'error')
            return False, 'error'
        snr_smooth = smooth_filter(self.snr, n_smooth)
        if snr_smooth is None:
            log('Could not smooth_filter self.snr for station "{}"'
                .format(self.station), 'error')
            log(self.snr, 'error')
            return False, 'error'
        self._set_snr_threshold(snr_smooth)
        sign_change = np.diff(np.sign(snr_smooth.data - self.snr_threshold))
        crossings = len(sign_change[sign_change == 2])
        if debug:
            snr_smooth.plot()
        max_cross = self.params.max_threshold_crossings
        trustworthy = (crossings > 0 and crossings <= max_cross)
        if not trustworthy:
            s = 'not trustworthy: '
            if crossings == 0:
                s += f'never crossed the threshold ({self.snr_threshold})'
            else:
                s += f'crossed the threshold ({self.snr_threshold}) '
                s += f'{crossings:d} times (> {max_cross:d})'
        else:
            s = "trustworthy"
        return (crossings > 0
                and crossings <= self.params.max_threshold_crossings), s

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
