# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# readmain.m


class PolarityParameters():
    """
    Polarity Parameters
    """
    def __init__(self,
                 calculate_window=2,
                 analyze_window=4,
                 DR_threshold_P=0.4,
                 DR_threshold_S=-0.4,
                 DR_smooth_length=1):
        """
        Initialize Polarity Parameters

        :param calculate_window: number of seconds AFTER a pick in which
            to calculate polarities
        :param analyze_window: number of seconds to analyze for each
            polarity calculation
        :param DR_threshold_P: DipRect P-wave threshold
        :param DR_threshold_S: DipRect S-wave threshold
        :param DR_smooth_length: seconds over which to smooth polarity and
            rectilinearity when calculating dip-rect
        """
        self.calculate_window = calculate_window
        self.analyze_window = analyze_window
        self.DR_threshold_P = DR_threshold_P
        self.DR_threshold_S = DR_threshold_S
        self.DR_smooth_length = DR_smooth_length

    def __str__(self):
        str = "PolarityParameters:\n"
        str += f"    calculate_window = {self.calculate_window}\n"
        str += f"    analyze_window = {self.analyze_window}\n"
        str += f"    DR_threshold_P = {self.DR_threshold_P}\n"
        str += f"    DR_threshold_S = {self.DR_threshold_S}\n"
        str += f"    DR_smooth_length = {self.DR_smooth_length}\n"
        return str

    @classmethod
    def from_dict(cls, thedict):
        return(cls(**thedict))


if __name__ == '__main__':
    pass
