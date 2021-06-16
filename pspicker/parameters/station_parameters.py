from .kurtosis_parameters import KurtosisParameters


class StationParameters():
    """
    Station Parameters
    """
    def __init__(self, picking_components, SNR_energy,
                 kurtosis, resp_file, use_polarity=False, lag=999,
                 max_candidates=5):
        """
        :param picking_components: PickingComponents dictionary
        :param SNR_energy: SNREnergy dictionary
        :param kurtosis: KurtosisParameters dictionary
        :param resp_file: response file name
        :param lag: unused
        :param use_polarity: Use polarity analysis to identify phases?
        :param max_candidates: maximum number of extrema to pick as candidates
        """
        self.picking_components = PickingComponents(**picking_components)
        self.SNR_energy = SNREnergy(**SNR_energy)
        self.kurtosis = KurtosisParameters(**kurtosis)
        self.resp_file = resp_file
        self.use_polarity = use_polarity
        self.lag = lag
        self.max_candidates = max_candidates

    @classmethod
    def from_dict(cls, thedict):
        return cls(**thedict)


class PickingComponents():
    """
    Station Parameter picking components
    """
    def __init__(self, P, S):
        """
        Components used to search for arrivals and amplitudes
        :param P: For P waves
        :param S: For S waves
        """
        self.P = P
        self.S = S

   #  @classmethod
   #  def from_dict(cls, thedict):
   #      return cls(**thedict)
      
  
class SNREnergy():
    """
    Station Parameter SNR and Energy parameters
    """
    def __init__(self, frequency_band, window):
        """
        :param frequency_band: frequency band in which to calculate
            SNR [lo, high]
        :param window: only look at data from t-winlen to t, where t
            is the time of the peak waveform energy.  If == 0, don't use
            energy criteria.
        """
        self.frequency_band = frequency_band
        self.window = window

    # @classmethod
    # def from_dict(cls, thedict):
    #     return cls(**thedict)


