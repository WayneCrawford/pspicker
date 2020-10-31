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
        assert P_write_cmp in "ZNEH", 'not in "ZNEH"'
        assert S_write_cmp in "ZNEH", 'not in "ZNEH"'
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
        """ Returns the full name of the component """
        return getattr(self, self._P_write_cmp)

    @property
    def S_write_cmp(self):
        """ Returns the full name of the component """
        return getattr(self, self._S_write_cmp)

    def __str__(self, format=None):
        """
        :param format: 'table_row': print as a table row
                       'table_header': print a table header (no object info)
                       None: print normally
        """
        if format == 'table_header':
            s = '{:^6s}|{:^6s}|{:^6s}|{:^6s}|{:^9s}|{:^9s}|{:^8s}|{:^8s}|'.\
                format('Z', 'N', 'E', 'H', 'P_write', 'S_write',
                       'P_phase', 'S_phase')
            # s += '{:-<6s}+{:-<6s}+{:-<6s}+{:-<6s}+'.format(
            #     '', '', '', '')
            # s += '{:-<9s}+{:-<9s}+{:-<8s}+{:-<8s}|'.format(
            #     '', '', '', '')

            return s
        elif format == 'table_row':
            if self.H is not None:
                H = self.H.split('.')[-1]
            else:
                H = 'None'
            s = '{:^6s}|{:^6s}|{:^6s}|{:^6s}|{:^9s}|{:^9s}|{:^8s}|{:^8s}|'.\
                format(self.Z.split('.')[-1], self.N.split('.')[-1],
                       self.E.split('.')[-1], H,
                       self._write_cmp_str(self._P_write_cmp),
                       self._write_cmp_str(self._S_write_cmp),
                       self.P_write_phase, self.S_write_phase)
            return s
        else:
            s = f'            Z: {self.Z}\n'
            s += f'            N: {self.N}\n'
            s += f'            E: {self.E}\n'
            s += f'            H: {str(self.H)}\n'
            s += f'  P_write_cmp: {self._P_write_cmp} (=>{self.P_write_cmp})\n'
            s += f'  S_write_cmp: {self._S_write_cmp} (=>{self.S_write_cmp})\n'
            s += f'P_write_phase: {self.P_write_phase}\n'
            s += f'S_write_phase: {self.S_write_phase}\n'
            return s

    def _write_cmp_str(self, write_cmp):
        return write_cmp + '(' + getattr(self, write_cmp).split('.')[-1] + ')'

    def _shortChName(self, write_to):
        """Shortened channel name of the given seed_id"""
        seed_id = self[write_to]     # Dict acces to class parameter
        return seed_id[-3] + seed_id[-1]
