class KurtosisParameters():
    """
    Kurtosis Parameters
    """
    def __init__(self, frequency_bands, window_lengths,
                 extrema_smoothings=[40], n_smooth=1):
        """
        :param frequency_bands: frequency bands to filter data before
            calculating kurtosis
        :type frequency_bands: list of [lo, high]s
        :param window_lengths:  window lengths (seconds) over which to
            calculate kurtosis
        :type window_lengths: list
        :param extrema_smoothings: kurtosis smoothing sequence (increasing
            list of samples)
        :type extrema_smoothings: list of ints
        :param n_smooth: samples over which to smooth kurtosis
        :type n_smooth: samples over which to smooth kurtosis
        """
        assert isinstance(frequency_bands, list)
        assert isinstance(window_lengths, list)
        assert isinstance(extrema_smoothings, list)
        for fb in frequency_bands:
            assert(len(fb) == 2)
        for x in window_lengths:
            assert isinstance(x, (int, float))
        for x in extrema_smoothings:
            assert isinstance(x, int)

        self.frequency_bands = frequency_bands
        self.window_lengths = window_lengths
        self.extrema_smoothings = extrema_smoothings
        self.n_smooth = n_smooth

    def __str__(self):
        str = "KurtosisParmeters("
        str += f"frequency_bands = {self.frequency_bands}, "
        str += f"window_lengths = {self.window_lengths}, "
        str += f"extrema_smoothings = {self.extrema_smoothings})"
        return str

    @classmethod
    def from_dict(cls, thedict):
        return cls(**thedict)
