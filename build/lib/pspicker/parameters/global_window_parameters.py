from .kurtosis_parameters import KurtosisParameters


class GlobalWindowParameters():
    """
    Global Window Parameters
    """
    def __init__(self,
                 kurtosis,
                 offsets,
                 distri_secs,
                 max_candidates=5,
                 end_cutoff=0.9):
        """
        Initialize Global Window Parameters

        :param kurtosis: KurtosisParameters dictionary
        :param distri_secs: Length in seconds of sliding window to use
            to find the densest pick time over all stations
        :param max_candidates: Number of pick candidates to use in evaluating
            the overall pick distribution
        :param offsets: Offsets in seconds from densest pick
            time for global (all station based) rewindowing [left, right]
        :param end_cutoff: What fraction of data (from start) to
            look at for global Kurtosis window. 1.0 looks everywhere
        """
        assert len(offsets) == 2, "len(offsets) != 2"

        self.kurtosis = KurtosisParameters(**kurtosis)
        assert len(self.kurtosis.frequency_bands) == 1
        self.offsets = offsets
        self.end_cutoff = end_cutoff
        self.distri_secs = distri_secs
        self.max_candidates = max_candidates

    def __str__(self):
        str = "GlobalWindowParameters:\n"
        str += f"    kurtosis = {self.kurtosis}\n"
        str += f"    distri_secs = {self.distri_secs}\n"
        str += f"    offsets = {self.offsets}\n"
        str += f"    end_cutoff = {self.end_cutoff}\n"
        str += f"    max_candidates = {self.max_candidates}\n"
        return str

    @classmethod
    def from_dict(cls, thedict):
        return cls(**thedict)


if __name__ == '__main__':
    pass
