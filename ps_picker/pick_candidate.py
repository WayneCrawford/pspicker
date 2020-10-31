"""
Pick candidate class
"""
from obspy.core import UTCDateTime
from obspy.core.event.origin import Pick, Arrival
from obspy.core.event.base import WaveformStreamID, QuantityError


class PickCandidate():
    def __init__(self, time, picker_type, picker_value, snr=None,
                 DR=None, phase_guess=None, sampling_rate=None, station=None,
                 weight=None):
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
        :param weight: Nordic-style pick weight [0=best, 3=worst]
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
        self.weight = weight

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

    def to_obspy(self, channel_maps, quality_thresholds=None):
        """
        Return obspy Pick and Arrival

        Arrival is None if there is no pick weight
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
        self.weight = self._get_weight(quality_thresholds)
        time_errors = self._time_error()
        pick = Pick(time=self.time,
                    time_errors=time_errors,
                    waveform_id=id,
                    phase_hint=phase_hint,
                    evaluation_mode='automatic',
                    evaluation_status='preliminary')
        if self.weight is not None:
            return pick, Arrival(pick_id=pick.resource_id,
                                 time_weight=self.weight)
        else:
            return pick, None

    def _get_weight(self, quality_thresholds):
        """
        Return pick weights corresponding to SNR

        :param quality_thresholds: list of 4 signal-to-noise ratio thresholds
            for pick quality, from lowest [quality=3] to highest [quality=0]
        :returns: pick quality (0 to 3)
        :rtype: obspy QuantityError
        """
        if quality_thresholds is None:
            return None
        quality_thresholds.sort()
        if self.snr > quality_thresholds[3]:
            return 0
        elif self.snr >= quality_thresholds[2]:
            return 1
        elif self.snr >= quality_thresholds[1]:
            return 2
        elif self.snr >= quality_thresholds[0]:
            return 3
        else:
            return 4

    def _time_error(self,):
        """
        Return approximate pick time errors corresponding to pick weight

        Errors are multiplied by 2 for self.phase_guess == 'S'
        :returns: time_errors
        :rtype: obspy QuantityError
        """
        if not isinstance(self.weight, int):
            return None
        elif self.sampling_rate is None:
            return None

        assert self.phase_guess in 'PS',\
            "phase_guess '{self.phase_guess}' not in 'PS'"

        if self.weight == 0:
            uncertainty = 2. / self.sampling_rate
        elif self.weight == 1:
            uncertainty = 8. / self.sampling_rate
        elif self.weight == 2:
            uncertainty = 32. / self.sampling_rate
        elif self.weight == 3:
            uncertainty = 128. / self.sampling_rate
        else:
            uncertainty = 2000. / self.sampling_rate
        if self.phase_guess == 'S':
            uncertainty *= 2.
        return QuantityError(uncertainty)
