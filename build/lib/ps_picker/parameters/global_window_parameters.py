# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# readmain.m


class GlobalWindowParameters():
    """
    Global Window Parameters
    """
    def __init__(self,
                 kurt_frequency_band,
                 kurt_window_length,
                 offsets,
                 distri_secs,
                 n_extrema=5,
                 kurt_extrema_smoothing=40,
                 end_cutoff=0.9):
        """
        Initialize Global Window Parameters

        :param kurt_frequency_band: The frequency band for kurtosis calculation
            during global rewindowing  [low, high]
        :param kurt_window_length: Length in seconds of sliding window for
            Kurtosis during global rewindowing
        :param distri_secs: Length in seconds of sliding window to use
            to find the densest pick time over all stations
        :param n_extrema: Number of "extrema" to use in evaluating the overall
            pick distribution
        :param kurt_extrema_smoothing: Number of samples to use in smoothing
            window when calculating extrema
        :param offsets: Offsets in seconds from densest pick
            time for global (all station based) rewindowing [left, right]
        :param end_cutoff: What fraction of data (from start) to
            look at for global Kurtosis window. 1.0 looks everywhere
        """
        assert len(kurt_frequency_band) == 2, "len(kurt_frequency_band) != 2"
        assert len(offsets) == 2, "len(offsets) != 2"

        self.kurt_frequency_band = kurt_frequency_band
        self.kurt_window_length = kurt_window_length
        self.offsets = offsets
        self.end_cutoff = end_cutoff
        self.distri_secs = distri_secs
        self.n_extrema = n_extrema
        self.kurt_extrema_smoothing = kurt_extrema_smoothing

    def __str__(self):
        str = "GlobalWindowParameters:\n"
        str += f"    kurt_frequency_band = {self.kurt_frequency_band}\n"
        str += f"    kurt_window_length = {self.kurt_window_length}\n"
        str += f"    distri_secs = {self.distri_secs}\n"
        str += f"    offsets = {self.offsets}\n"
        str += f"    end_cutoff = {self.end_cutoff}\n"
        str += f"    n_extrema = {self.n_extrema}\n"
        str += f"    kurt_extrema_smoothing = {self.kurt_extrema_smoothing}\n"
        return str

    @classmethod
    def from_dict(cls, thedict):
        return cls(**thedict)


if __name__ == '__main__':
    pass
