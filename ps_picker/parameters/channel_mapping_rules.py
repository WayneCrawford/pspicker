# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# readmain.m


class ChannelMappingRules:
    """
    Rules for mapping channels
    """
    def __init__(self, compZ='Z3', compE='E2X', compN='N1Y', compH='HF',
                 S_write_cmp='N', P_write_cmp='Z',
                 P_write_phase='Pg', S_write_phase='Sg',
                 band_order='GFDCEHSBMLV'):
        """
        :param compZ: comp characters to assign to Z component
        :param compE: comp characters to assign to E component
        :param compN: comp characters to assign to N component
        :param compH: comp characters to assign to H component
        :param S_write_cmp: save S picks to this component.  Should be
            one of 'ZNEH': output will be written the the component matching
            this value using the comp? regexp above
        :param P_write_cmp: save P picks to this component.  Should be
            one of 'ZNEH'
        :param S_write_phase: call S picks this phase name
        :param P_write_phase: call P picks this phase name
        :param band_order: order of preference for band codes (if component
            choice does not yield unique selection)
        """
        assert S_write_cmp in 'ZNEH',\
            f"S_write_cmp = {S_write_cmp}, not in 'ZNEH'"
        assert P_write_cmp in 'ZNEH',\
            f"P_write_cmp = {P_write_cmp}, not in 'ZNEH'"
        self.compZ = compZ
        self.compE = compE
        self.compN = compN
        self.compH = compH
        self.S_write_cmp = S_write_cmp
        self.P_write_cmp = P_write_cmp
        self.P_write_phase = P_write_phase
        self.S_write_phase = S_write_phase
        self.band_order = band_order

    def __str__(self):
        str = "ChannelMappingRules\n"
        str += f"    compZ = {self.compZ}\n"
        str += f"    compE = {self.compE}\n"
        str += f"    compN = {self.compN}\n"
        str += f"    compH = {self.compH}\n"
        str += f"    S_write_cmp = {self.S_write_cmp}\n"
        str += f"    P_write_cmp = {self.P_write_cmp}\n"
        str += f"    P_write_phase = {self.P_write_phase}\n"
        str += f"    S_write_phase = {self.S_write_phase}\n"
        return str

    @classmethod
    def from_dict(cls, thedict):
        return(cls(**thedict))
