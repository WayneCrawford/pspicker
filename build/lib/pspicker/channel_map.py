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
        assert Z is not None
        self.Z = Z
        self.N = N
        self.E = E
        self.H = H
        self._P_write_cmp = P_write_cmp
        self._S_write_cmp = S_write_cmp
        self.P_write_phase = P_write_phase
        self.S_write_phase = S_write_phase

    @property
    def P_write_seed_string(self):
        """ Returns the full seed string name"""
        return getattr(self, self._P_write_cmp)

    @property
    def S_write_seed_string(self):
        """ Returns the full seed string """
        return getattr(self, self._S_write_cmp)

    @property
    def P_write_cmp(self):
        """ Returns the component name"""
        return self.P_write_seed_string.split('.')[-1]

    @property
    def S_write_cmp(self):
        """ Returns the component name"""
        return self.S_write_seed_string.split('.')[-1]

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
            H, Z, N, E = 'None', 'None', 'None', 'None'
            if self.H is not None:
                H = self.H.split('.')[-1]
            if self.Z is not None:
                Z = self.Z.split('.')[-1]
            if self.N is not None:
                N = self.N.split('.')[-1]
            if self.E is not None:
                E = self.E.split('.')[-1]
            s = '{:^6s}|{:^6s}|{:^6s}|{:^6s}|{:^9s}|{:^9s}|{:^8s}|{:^8s}|'.\
                format(Z, N, E, H,
                       self.P_write_cmp, self.S_write_cmp,
                       self.P_write_phase, self.S_write_phase)
#                       self._write_cmp_str(self._P_write_cmp),
#                       self._write_cmp_str(self._S_write_cmp),
#                       self.P_write_phase, self.S_write_phase)
            return s
        else:
            s = f'            Z: {str(self.Z)}\n'
            s += f'            N: {str(self.N)}\n'
            s += f'            E: {str(self.E)}\n'
            s += f'            H: {str(self.H)}\n'
            s += f'  P_write_cmp: {self._P_write_cmp} (=>{self.P_write_cmp})\n'
            s += f'  S_write_cmp: {self._S_write_cmp} (=>{self.S_write_cmp})\n'
            s += f'P_write_phase: {self.P_write_phase}\n'
            s += f'S_write_phase: {self.S_write_phase}\n'
            return s

    def _shortChName(self, write_to):
        """Shortened channel name of the given seed_id"""
        seed_id = self[write_to]     # Dict acces to class parameter
        return seed_id[-3] + seed_id[-1]
