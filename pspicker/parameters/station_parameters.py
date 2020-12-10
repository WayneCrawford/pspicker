from .kurtosis_parameters import KurtosisParameters


class StationParameters():
    """
    Station Parameters
    """
    def __init__(self, P_comp, S_comp, energy_frequency_band, energy_window,
                 kurtosis, resp_file, use_polarity=False, lag=999,
                 max_candidates=5):
        """
        :param P_Comp: String of components used to search for P arrival
            & amplitude
        :param S_comp: String of components used to search for S arrival
            & amplitude
        :param energy_frequency_band: frequency band in which to calculate
            SNR [lo, high]
        :param energy_window: only look at data from t-winlen to t, where t
            is the time of the peak waveform energy.  If == 0, don't use
            energy criteria.
        :param kurtosis: KurtosisParameters dictionary
        :param resp_file: response file name
        :param lag: unused
        :param use_polarity: Use polarity analysis to identify phases?
        :param max_candidates: maximum number of extrema to pick as candidates
        """
        self.P_comp = P_comp
        self.S_comp = S_comp
        self.energy_frequency_band = energy_frequency_band
        self.kurtosis = KurtosisParameters(**kurtosis)
        self.energy_window = energy_window
        self.resp_file = resp_file
        self.use_polarity = use_polarity
        self.lag = lag
        self.max_candidates = max_candidates

    @classmethod
    def from_dict(cls, thedict):
        return cls(**thedict)
