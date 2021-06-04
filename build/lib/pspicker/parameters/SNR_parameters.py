# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# readmain.m


class SNRParameters():
    """
    Picker Parameters
    """
    def __init__(self,
                 signal_window,
                 noise_window,
                 quality_thresholds,
                 max_threshold_crossings=5,
                 threshold_parameter=0.2):
        """
        Initialize signal-to-noise ratio parameters

        :param noise_window: Length of noise window in seconds
        :param signal_window: Length of signal window in seconds
        :param quality_thresholds: SNR thresholds for pick quality
            list with values for quality = '3', '2', '1' and '0', in order
            the minimum value is also used as a baseline for SNR_threshold
        :param threshold_parameter: Threshold for SNR-based quality
            evaluation.
            if > 0 and less than 1, then SNR_threshold is this times the
                maximum SNR
            if < 0, then SNR_threshold is abs(SNR_threshold_parameter)
        :param max_threshold_crossings: Maximum crossings of SNR threshold
            level to accept a trace
        """
        qt = quality_thresholds
        assert len(qt) == 4, f'Need 4 quality_thresholds, found {len(qt):d}'
        assert qt == sorted(qt), f"quality thresholds not increasing: {qt}"

        self.signal_window = signal_window
        self.noise_window = noise_window
        self.quality_thresholds = quality_thresholds
        self.max_threshold_crossings = max_threshold_crossings
        self.threshold_parameter = threshold_parameter

    def __str__(self):
        str = "SNRParameters:\n"
        str += f"    signal_window = {self.signal_window}\n"
        str += f"    noise_window = {self.noise_window}\n"
        str += f"    quality_thresholds = {self.quality_thresholds}\n"
        str += "    max_threshold_crossings = "
        str += f"{self.max_threshold_crossings}\n"
        str += f"    threshold_parameter = {self.threshold_parameter}\n"
        return str

    @classmethod
    def from_dict(cls, thedict):
        return cls(**thedict)


if __name__ == '__main__':
    pass
