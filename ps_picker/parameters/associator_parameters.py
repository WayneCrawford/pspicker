# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# readmain.m


class AssociatorParameters():
    """
    Associator Parameters
    """
    def __init__(self,
                 cluster_window_P,
                 cluster_window_S,
                 distri_min_values=4,
                 distri_nstd_picks=3.2,
                 distri_nstd_delays=4):
        """
        Initialize Associator Parameters

        :param cluster_window_P: Window [s] for P-pick rejection
                                       based on clustering
        :param cluster_window_S: Same as above, for S-pick
        :param distri_min_values: Minimum number of values (P-picks,
            S-picks or PS delays) needed to perform distribution-based
            rejection (if > n_stations of n_picks, will not perform
            distribution-based rejection
        :param distri_nstd_picks: maximum number of deviations from
            standard distribution to accept for P picks, and for S picks
        :param distri_nstd_delays: same as above, for P-S delays
        """
        self.cluster_window_P = cluster_window_P
        self.cluster_window_S = cluster_window_S
        self.distri_min_values = float(distri_min_values)
        self.distri_nstd_picks = distri_nstd_picks
        self.distri_nstd_delays = distri_nstd_delays

    def __str__(self):
        str = "AssociatorParameters:\n"
        str += f"    cluster_window_P = {self.assoc_cluster_window_P}\n"
        str += f"    cluster_window_S = {self.assoc_cluster_window_S}\n"
        str += f"    distri_min_values = {self.distri_min_values}\n"
        str += f"    distri_nstd_picks = {self.distri_nstd_picks}\n"
        str += f"    distri_nstd_delays = {self.distri_nstd_delays}\n"
        return str

    @classmethod
    def from_dict(cls, thedict):
        return cls(**thedict)


if __name__ == '__main__':
    pass
