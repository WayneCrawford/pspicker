# Generated with SMOP  0.41-beta
# from smop.libsmop import *


class StationParameters():
    """
    Station Parameters
    """
    def __init__(self, P_comp, S_comp, energy_frequency_band,
                 kurt_frequency_bands, kurt_window_lengths,
                 kurt_extrema_smoothings, energy_window, resp_file,
                 use_polarity=False, lag=999, n_extrema=5):
        """
        param P_Comp: String of components used to search for P arrival
            & amplitude
        param S_comp: String of components used to search for S arrival
            & amplitude
        param energy_frequency_band: frequency band in which to calculate
            SNR [lo, high]
        :param kurt_frequencies: kurtosis frequency bands (list of [lo, high]s)
        :param kurt_window_lengths: kurtosis window lengths (list of seconds)
        :param kurt_extrema_smoothings: kurtosis smoothing sequence (increasing
            list of samples)
        :param energy_window: only look at data from t-winlen to t, where t
            is the time of the peak waveform energy.  If == 0, don't use
            energy criteria.
        :param resp_file: response file name
        :param lag: unused
        :param use_polarity: Use polarity analysis to identify phases?
        :param n_extrema: maximum number of extrema to pick as candidates
        """
        self.P_comp = P_comp
        self.S_comp = S_comp
        self.energy_frequency_band = energy_frequency_band
        self.kurt_frequency_bands = kurt_frequency_bands
        self.kurt_window_lengths = kurt_window_lengths
        self.kurt_extrema_smoothings = kurt_extrema_smoothings
        self.energy_window = energy_window
        self.resp_file = resp_file
        self.use_polarity = use_polarity
        self.lag = lag
        self.n_extrema = n_extrema

    @classmethod
    def from_dict(cls, thedict):
        return cls(**thedict)
