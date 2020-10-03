# from readMSEEDTraces.m
class ChannelMap():
    """
    Class to store seed_ids corresponing to diferent components
    """
    def __init__(self, Z, N, E, H, P_write_cmp, P_write_phase,
                 S_write_cmp, S_write_phase):
        """
        :paran Z: seed_id corresponding to station's vertical component
        :paran N: seed_id corresponding to station's north component
        :paran E: seed_id corresponding to station's east component
        :paran H: seed_id corresponding to station's hydrophone/pressure
            component
        :param P_write_cmp: component to write P picks to
        :param S_write_cmp: component to write S picks to
        :param P_write_phase: what phase to name P picks
        :param S_write_phase: what phase to name S picks
        """
        self.Z = Z
        self.N = N
        self.E = E
        self.H = H
        self._P_write_cmp = P_write_cmp
        self._S_write_cmp = S_write_cmp
        self.P_write_phase = P_write_phase
        self.S_write_phase = S_write_phase

    @property
    def P_write_cmp(self):
        return getattr(self, self._P_write_cmp)

    @property
    def S_write_cmp(self):
        return getattr(self, self._S_write_cmp)

    def __str__(self):
        fmt = '{:10s} | Z: {:-10s} | N: {:-10s} | E: {:-10s} | H: {:-10s} |'
        return fmt.format('', self.Z, self.N, self.E, self.H)

    def _shortChName(self, write_to):
        """Shortened channel name of the given seed_id"""
        seed_id = self[write_to]     # Dict acces to class parameter
        return seed_id[-3] + seed_id[-1]
