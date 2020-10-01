# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# readmain.m


class ChannelMappingRules:
    """
    Rules for mapping channels
    """
    def __init__(self, compZ='[SHDE][HL][Z3]', compE='[SHDE][HL][XE2]',
                 compN='[SHDE][HL][YN1]', compH='([SHDE]DH|P|SDF)',
                 putPick_S_Comp='N', putPick_P_Comp='Z',
                 putPick_P_Phase='Pg', putPick_S_Phase='Sg'):
        """
        :param compZ: channel name regexp to assign to Z component
        :param compE: channel name regexp to assign to E component
        :param compN: channel name regexp to assign to N component
        :param compH: channel name regexp to assign to H component
        :param putPick_S_Comp: save S picks to this component
        :param putPick_P_Comp: save P picks to this component
        :param putPick_S_Phase: call S picks this phase name
        :param putPick_P_Phase: call P picks this phase name

        The channel name regexps must match the end of the component name.
        For example:
            '[XE2]'         accepts names 'HDX', 'BHX', 'BHE', 'BH2'...
            '(X|E|2)'       does the same
            'SHX'           accepts the name 'SHX'
            '(SHX|SHE|SH2)'	accepts the names 'SHX', 'SHE' or 'SH2'
            ''              accepts nothing
        """
        self.compZ = compZ
        self.compE = compE
        self.compN = compN
        self.compH = compH
        self.putPick_S_Comp = putPick_S_Comp
        self.putPick_P_Comp = putPick_P_Comp
        self.putPick_P_Phase = putPick_P_Phase
        self.putPick_S_Phase = putPick_S_Phase

    def __str__(self):
        str = "ChannelMappingRules\n"
        str += f"    compZ = {self.compZ}\n"
        str += f"    compE = {self.compE}\n"
        str += f"    compN = {self.compN}\n"
        str += f"    compH = {self.compH}\n"
        str += f"    putPick_S_Comp = {self.putPick_S_Comp}\n"
        str += f"    putPick_P_Comp = {self.putPick_P_Comp}\n"
        str += f"    putPick_P_Phase = {self.putPick_P_Phase}\n"
        str += f"    putPick_S_Phase = {self.putPick_S_Phase}\n"
        return str
