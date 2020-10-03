# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# readmain.m


# class OtherParameters():
#    """
#    Other parameters that were set for pickfunction
#
#    Once I figure out what/why they are, give this class a name and/or
#    subdivide it
#    """
#    def __init__(self, global_pick_wind_end=90, cleandistri_nPicks=1.e9,
#                 cleandistri_PorS=3.2, cleandistri_PtoS=4,
#                 responseFileType=''):
#        """
#        :param global_pick_wind_end: What percentage of data (from start) to
#            look at for global Kurtosis window. 100 looks everywhere
#        :param cleandistri_nPicks: Minimum number of picks needed for
#            distribution-based rejection
#        :param cleandistri_PorS: maximum number of deviations from standard
#            distribution to accept
#        :param cleandistri_PtoS: ????
#        :param responseFileType: Format of the response file(s).  'GSE' or
#            empty.  If empty use weird Baillardesque PoleZero format.
#        """
#        self.global_pick_wind_end = global_pick_wind_end
#        self.cleandistri_nPicks = cleandistri_nPicks
#        self.cleandistri_PorS = cleandistri_PorS
#        self.cleandistri_PtoS = cleandistri_PtoS
#        self.responseFileType = responseFileType
#
#    def __str__(self):
#        str = "OtherParameters\n"
#        str += f"    global_pick_wind_end = {self.global_pick_wind_end}\n"
#        str += f"    cleandistri_nPicks = {self.cleandistri_nPicks}\n"
#        str += f"    cleandistri_PorS = {self.cleandistri_PorS}\n"
#        str += f"    cleandistri_PtoS = {self.cleandistri_PtoS}\n"
#        str += f"    responseFileType = {self.responseFileType}\n"
#        return str
