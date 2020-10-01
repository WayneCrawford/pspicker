# Generated with SMOP  0.41-beta
# from smop.libsmop import *


class StationParameters():
    """
    Station Parameters
    """
    def __init__(self, P_comp, S_comp, f_energy, frequencies, window_lengths,
                 smoothing_sequence, energy_window, response_file,
                 use_polarity=False, lag=999, n_follow=2):
        """
        param P_Comp: String of components used to search for P arrival
            & amplitude
        param S_comp: String of components used to search for S arrival
            & amplitude
        param f_energy:  frequency range in which to search for SNR [lo, high]
        :param frequencies: kurtosis frequency bands (list of [lo, high]s)
        :param window_lengths: kutosis window lengths (list of seconds)
        :param smoothing_sequence: kurtosis smooting sequence (increasing
            list of samples(?))
        :param energy_window: only look at data from t-winlen to t, where t
            is the time of the peak waveform energy.  If == 0, don't use
            energy criteria.
        :param response: response file name
        :param lag: unused
        :param use_polarity: Use polarity analysis to identify phases?
        :param n_follow: number of extrema to follow (1 or 2).  Generally use
            2 (S and P) unless data are problematic (some OBS data,
            for example)
        """
        self.P_comp = P_comp
        self.S_comp = S_comp
        self.f_energy = f_energy
        self.frequencies = frequencies
        self.window_lengths = window_lengths
        self.smoothing_sequence = smoothing_sequence
        self.energy_window = energy_window
        self.response_file = response_file
        self.use_polarity = use_polarity
        self.lag = lag
        self.n_follow = n_follow
