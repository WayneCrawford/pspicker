"""
Pick candidate class
"""
from obspy.core import UTCDateTime
from obspy.core.event.origin import Pick
from obspy.core.event.base import WaveformStreamID, QuantityError


class PickCandidate():
    def __init__(self, time, picker_type, picker_value, snr=None,
                 DR=None, phase_guess=None, sampling_rate=None, station=None):
        """
        :param time: time of the PickCandidate
        :param picker_type: picker type (e.g. 'kurtosis')
        :param picker_value: measure of pick quality (higer=better) made
            by the picker
        :param snr: signal-to-noise ratio at the pick
        :param DR: Dip-Rect parameter associated with the PickCandidate
        :param phase_guess: most likely phase type
        :param station: station for this candidate
        :param sampling_rate: sampling rate of data used to make the candidate
        """
        assert isinstance(time, UTCDateTime), 'time is not a UTCDateTime'
        assert isinstance(picker_type, str), 'picker_type is not a str'
        picker_value = float(picker_value)

        self.time = time
        self.picker_type = picker_type
        self.picker_value = picker_value
        self.snr = snr
        self.DR = DR
        self.phase_guess = phase_guess
        self.sampling_rate = sampling_rate
        self.station = station

    def __str__(self):
        s = f"PickCandidate: {self.time.strftime('%Y%m%dT%H%M%S')}"
        s += f", {self.picker_type} value={self.picker_value:.3g}"
        if self.snr is None:
            s += ", snr=None"
        else:
            s += f", snr={self.snr:.3g}"
        if self.DR is not None:
            s += f", DR={self.DR:.3g}"
        if self.phase_guess is not None:
            s += f', phase_guess="{self.phase_guess}"'
        if self.sampling_rate is not None:
            s += f', sampling_rate={self.sampling_rate:.3g}'
        if self.station is not None:
            s += f', station="{self.station}"'
        return s

    def to_pick(self, channel_maps, quality_thresholds=None):
        """
        Return obspy Pick
        :param channel_maps: {station: ChannelMap object}
        :param quality_thresholds: list of 4 signal-to-noise ratio thresholds
            for pick quality, from lowest [quality=3] to highest [quality=0]
        :parm sampling_rate: sampling rate (used to estimate timing error)
        """
        c_map = channel_maps[self.station]
        if self.phase_guess == 'P':
            id = WaveformStreamID(seed_string=c_map.P_write_cmp)
            phase_hint = c_map.P_write_phase
        elif self.phase_guess == 'S':
            id = WaveformStreamID(seed_string=c_map.S_write_cmp)
            phase_hint = c_map.S_write_phase
        else:
            raise ValueError("phase guess '{self.phase_guess}' not 'P' or 'S'")
        if quality_thresholds:
            time_errors = self._time_error(quality_thresholds)
        else:
            time_errors = None
        return Pick(time=self.time,
                    time_errors=time_errors,
                    waveform_id=id,
                    phase_hint=phase_hint,
                    evaluation_mode='automatic',
                    evaluation_status='preliminary')

    def _time_error(self, quality_thresholds):
        """
        Return approximate pick time errors corresponding to SNR

        Errors are multiplied by 2 for self.phase_guess == 'S'
        :param quality_thresholds: list of 4 signal-to-noise ratio thresholds
            for pick quality, from lowest [quality=3] to highest [quality=0]
        :returns: time_errors
        :rtype: obspy QuantityError
        """
        if self.sampling_rate is None:
            return None
        assert self.phase_guess in 'PS',\
            "phase_guess '{self.phase_guess}' not in 'PS'"
        if self.snr > quality_thresholds[3]:
            uncertainty = 2. / self.sampling_rate
        elif self.snr >= quality_thresholds[2]:
            uncertainty = 8. / self.sampling_rate
        elif self.snr >= quality_thresholds[1]:
            uncertainty = 32. / self.sampling_rate
        elif self.snr >= quality_thresholds[0]:
            uncertainty = 128. / self.sampling_rate
        else:
            uncertainty = 2000. / self.sampling_rate
        if self.phase_guess == 'S':
            uncertainty *= 2.
        return QuantityError(uncertainty)
