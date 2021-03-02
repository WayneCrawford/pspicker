# import warnings
import itertools

import numpy as np
from scipy.cluster.hierarchy import fclusterdata
from obspy.taup import TauPyModel
from obspy.taup.taup_create import build_taup_model
from obspy.taup.tau_model import TauModel
from obspy.core import UTCDateTime
# from obspy.core.event.origin import Pick

from .utils import picks_matched_stations
from .logger import log
from .pick_candidate import PickCandidate


class Associator():
    """
    Associate Picks

    Currently only rejects picks that deviate from expected distribution
    ToDo:
        - Replace rejected picks with other candidates
        - Add velocity model-based selection criteria
    """
    # Values from Kennett 1995 for 0, 1, 5, 10, 30 and 90 degrees
    # default_delays = {'op': [0, 19.17, 76.27, 144.90, 370.27, 781.35],
    #                   'ps': [0, 13.92, 59.63, 114.2, 300.00, 654.43]}

    def __init__(self, params, vel_model=None, target_depth=10., max_offset=1.,
                 verbose=True, debug=False):
        """
        :param params: AssociatorParameters object
        :param vel_model: velocity model (obspy taup valid string or the
            absolute path name of a .nd file)
        :param target_depth: event depth to use to construct travel-time tables
        :param max_offset: maximum offset(degrees) for which to calculate
            travel_times

        The nd file format allows c-style comment lines, depth-parameter lines
        and named discontinuity lines.
        The depth-parameter lines have space separated columns for:
        Depth(km), Pvel(km/s), Svel(km/s) and, optional density(g/cc), Qp & Qs
        The named discontinuity lines are placed between the depth-parameter
        lines corresponding to the values above and below.  For example:
        // This is a comment
        0.0    5.0   3.0
        20     5.0   3.0
        20     6.5   3.7
        33     6.5   3.7
        mantle        # the word "mantle" designates that this is the moho
        33     7.8   4.4
        410    8.9   4.7
        """
        self.method = params.method
        self.distri_min_values = params.distri_min_values
        self.cluster_window_P = params.cluster_window_P
        self.cluster_window_S = params.cluster_window_S
        self.distri_nstd_picks = params.distri_nstd_picks
        self.distri_nstd_delays = params.distri_nstd_delays
        self.cluster_window_otime = params.cluster_window_otime
        self.vp_over_vs = params.otime_vp_vs
        self.p_cluster = {}   # a dictionary with key=station and value=ptime
        self.s_cluster = {}   # a dictionary with key=station and value=stime
        self.ps_cluster = {}   # a dictionary with key=station and value=stime
        self.o_cluster = {}   # a dictionary with key=station and value=otime
        self.delays = None
        if vel_model is not None:
            try:
                self.delays = self._get_delays(vel_model, target_depth,
                                               max_offset)
            except Exception:
                log("Couldn't generate ps delays using '" + vel_model
                    + "', using default", "warning")
        self.verbose = verbose
        self.debug = debug
        self.min_clust = 3  # minimum number of values in a "cluster"

    def run(self, picks, candidates):
        """
        Run associator

        Tries to associate by origin time, if that doesn't work associates by
        pick clustering
        :param picks: list of preferred picks
        :param candidates: all pick candidates (including preferred picks)

        Future upgrade: associate by position (requires velocity model and
        station positions)
        """
        picks = PickCandidate.remove_duplicates(picks)
        input_picks = picks.copy()
        method = self.method
        assert method in ['origin_time', 'arrival_time']
        if method == 'origin_time':
            assoc_str = 'origin times'
            picks, associated = self.find_same_origin_time(picks, candidates)
            if not associated:
                log('Could not associate by origin time', 'verbose')
                method = 'pick_time'
        if method == 'arrival_time':
            assoc_str = 'pick clustering'
            picks = self.remove_nonclustered(picks)
        picks = PickCandidate.remove_duplicates(picks)
        self._assoc_stats(assoc_str, input_picks, picks)
        return picks

    @staticmethod
    def _assoc_stats(assoc_type, input_picks, output_picks):
        """
        Print information about accepted, rejected, added and changed picks
        """
        ops = output_picks.copy()
        stats = {'accepted': [], 'rejected': [], 'modified': [], 'added': []}
        # n_changes = 0
        for p in input_picks:
            if p in ops:
                stats['accepted'].append(p.shortname)
                ops.remove(p)
            else:
                op_remove = None
                for op in ops:
                    if (p.station == op.station
                            and p.phase_guess == op.phase_guess):
                        stats['modified'].append(p.shortname)
                        op_remove = op
                if op_remove is not None:
                    ops.remove(op_remove)
                else:
                    stats['rejected'].append(p.shortname)
        for p in ops:
            stats['added'].append(p.shortname)

        s = f'Associated by {assoc_type}: '
        for key in stats:
            s += f'{len(stats[key]):d} {key}, '
        log(s, 'verbose')
        for key in ['rejected', 'modified', 'added']:
            if len(stats[key]) > 1:
                log(f"{key.upper()}: {stats[key]}", 'debug')

    def remove_nonclustered(self, picks):
        """
        Very basic (clustering-based) pick selection by association

        :param picks: list of obspy Picks
        :candidates: dictionary with key=station, value=list of PickCandidates
        """
        # Remove lines with arrival times = start_time
        # P_picked_cell=rm_cell_line(P_picked_cell[2] == 0,P_picked_cell)
        # S_picked_cell=rm_cell_line(S_picked_cell[2] == 0,S_picked_cell)

        picks = self._remove_unclustered(picks)
        stations = _pick_stations(picks)
        if self.distri_min_values <= len(stations):
            picks = self._remove_badly_distributed(picks)
            picks = self._remove_bad_delays(picks)
        return picks

    def find_same_origin_time(self, picks, candidates):
        """
        Select picks by origin time

        :param picks: list of PickCandidates pre-selected as P and S picks
        :candidates: list of all PickCandidates

        :returns: list of picks, whether associator worked
        :rtype: list of PickCandidate, bool
        """
        stations = _pick_stations(picks)
        ots = []
        for station in stations:
            station_picks = [x for x in picks
                             if x.station == station]
            p_pick = [x for x in station_picks if x.phase_guess == "P"]
            s_pick = [x for x in station_picks if x.phase_guess == "S"]
            if not (len(p_pick) == 1 and len(s_pick) == 1):
                continue
            ots.append(self._ps_to_otime(p_pick[0].time, s_pick[0].time))
        if len(ots) < self.min_clust:
            log(f'less than {self.min_clust} P-S calculated origin times'
                ', cannot associate by this criteria', 'verbose')
            for p in picks:
                log(f' {p}', 'verbose')
            return picks, False
        good_ots = self._cluster_clean_otimes([t for t in ots])
        if len(good_ots) < self.min_clust:
            log('less than {} P-S origin times agree, cannot associate'
                .format(self.min_clust), 'verbose')
            for p in picks:
                log(f' {p}', 'debug')
            return picks, False
        mean_ot = UTCDateTime(np.mean([x.timestamp for x in good_ots]))
        new_picks = self._find_otime_matching(mean_ot, picks, candidates)
        return new_picks, True

    def _find_otime_matching(self, ot, picks, candidates):
        """
        Return picks or candidates whose P-S delay matches the origin time

        :param ot: desired origin time
        :param picks: list of preferred PickCandidates
        :param candidates: list of all PickCandidates
        :returns: new list of preferred PickCandidates
        """
        otime_margin = self.cluster_window_otime / 2.
        stations_c = _pick_stations(candidates)
        # Error check
        for sta in _pick_stations(picks):
            if sta not in stations_c:
                raise ValueError(f"picked station {sta} not in candidates")

        new_picks = []
        self.o_cluster = {}
        for station in stations_c:
            sta_candidates = [x for x in candidates if x.station == station]
            p_existing, s_existing = self._get_existing(picks, station)

            # If we have a P and an S candidate that match ot, keep them
            if s_existing is not None and p_existing is not None:
                # log(f'{station}: s_existing and p_existing', 'debug')
                otime = self._ps_to_otime(p_existing.time, s_existing.time)
                if np.abs(otime - ot) < otime_margin:
                    new_picks.extend([p_existing, s_existing])
                    self.o_cluster[station] = otime
                    continue
            # Otherwise, if a P cand matches another cand's ot, keep them
            if p_existing is not None:
                # log(f'{station}: p_existing top', 'debug')
                # log(p_existing, 'debug')
                cands = [x for x in sta_candidates if not x == p_existing]
                if len(cands) > 0:
                    otimes = [self._ps_to_otime(p_existing.time, x.time)
                              for x in cands]
                    offsets = [np.abs(x - ot) for x in otimes]
                    if np.min(offsets) < otime_margin:
                        s_candidate = cands[np.argmin(offsets)]
                        s_candidate.phase_guess = 'S'
                        # log(p_existing, 'debug')
                        # log(s_candidate, 'debug')
                        new_picks.extend([p_existing, s_candidate])
                        self.o_cluster[station] = otimes[np.argmin(offsets)]
                        continue
            # Otherwise, an S cand that matches with another cand, keep them
            if s_existing is not None:
                # log(f'{station}: s_existing top', 'debug')
                cands = [x for x in sta_candidates if not x == s_existing]
                if len(cands) > 0:
                    otimes = [self._ps_to_otime(x.time, s_existing.time)
                              for x in cands]
                    offsets = [np.abs(x - ot) for x in otimes]
                    if np.min(offsets) < otime_margin:
                        p_candidate = cands[np.argmin(offsets)]
                        p_candidate.phase_guess = 'P'
                        new_picks.extend([p_candidate, s_existing])
                        self.o_cluster[station] = otimes[np.argmin(offsets)]
                        continue
            # Otherwise, if any combination of cands matches ot, keep them
            if sta_candidates is not None:
                # log(f'{station}: neither s_existing nor p_existing', 'debug')
                sort_cands = sorted(sta_candidates, key=lambda x: x.time)
                min_offset = otime_margin
                new_p, new_s = None, None
                for p_maybe, s_maybe in itertools.combinations(sort_cands, 2):
                    otime = self._ps_to_otime(p_maybe.time, s_maybe.time)
                    offset = np.abs(otime - ot)
                    if offset < min_offset:
                        new_p = p_maybe
                        new_s = s_maybe
                        new_otime = otime
                        min_offset = offset
                if new_p is not None:
                    new_p.phase_guess = 'P'
                    new_s.phase_guess = 'S'
                    new_picks.extend([new_p, new_s])
                    self.o_cluster[station] = new_otime
                    continue
            # Otherwise, keep a solitary P or S pick
            # (if I used station positions, I could compare times here)
            if p_existing is not None and s_existing is None:
                # log(f'{station}: p_existing bottom', 'debug')
                assert isinstance(p_existing, PickCandidate)
                new_picks.append(p_existing)
            elif p_existing is None and s_existing is not None:
                # log(f'{station}: s_existing bottom', 'debug')
                assert isinstance(s_existing, PickCandidate)
                new_picks.append(s_existing)
        # log('about to test for duplicates', 'debug')
        # new_picks = PickCandidate.remove_duplicates(new_picks)
        # log('just tested for duplicates', 'debug')
        return new_picks

    @staticmethod
    def _get_existing(picks, station):
        p_existings = [x for x in picks if x.station == station
                       and x.phase_guess == 'P']
        s_existings = [x for x in picks if x.station == station
                       and x.phase_guess == 'S']
        p_existing, s_existing = None, None
        if len(p_existings) > 0:
            p_existing = p_existings[0]
        if len(s_existings) > 0:
            s_existing = s_existings[0]
        return p_existing, s_existing

    def _remove_unclustered(self, picks):
        """
        Remove picks and reassign phases based on clustering

        Very basic association
        :param picks: list of PickCandidate objects
        :returns: modified list
        """
        # Pick_Function.m:746
        p_picks = [p for p in picks if p.phase_guess == 'P']
        if len(p_picks) >= self.distri_min_values:
            p_picks = self._cluster_clean_picks('P'. p_picks)
            log(f'clustered p_picks = {p_picks}', 'debug')
            self.p_cluster = {x.station: x.time for x in p_picks}

        s_picks = [p for p in picks if p.phase_guess == 'S']
        if len(s_picks) >= self.distri_min_values:
            s_picks = self._cluster_clean_picks('S', s_picks)
            log(f'clustered s_picks = {s_picks}', 'debug')
            self.p_cluster = {x.station: x.time for x in s_picks}

        return p_picks + s_picks

    def _remove_badly_distributed(self, picks):
        """
        Remove picks that are well outside of pick distribution

        :param picks: input Picks
        """
        p_picks = [p for p in picks if p.phase_guess == 'P']
        s_picks = [p for p in picks if p.phase_guess == 'S']
        _, iP = clean_distri([x.time.timestamp for x in p_picks],
                             self.distri_nstd_picks, 'median',
                             self.distri_min_values)
        _, iS = clean_distri([x.time.timestamp for x in s_picks],
                             self.distri_nstd_picks, 'median',
                             self.distri_min_values)
        return [p_picks[i] for i in iP] + [s_picks[i] for i in iS]

    def _remove_bad_delays(self, picks):
        """
        Eliminate picks based on distribution of P-S delays
        """
        # Pick_Function.m:763
        matches = picks_matched_stations(picks)
        if len(matches) == 0:
            return picks
        delay_stations = [m['station'] for m in matches]
        delays = [m['pickS'].time - m['pickP'].time for m in matches]
        _, i_PS = clean_distri(delays, self.distri_nstd_delays,
                               'median', self.distri_min_values)
        # Include good delay stations AND non-delay stations
        good_delay_stations = [delay_stations[i] for i in i_PS]
        good_delays = [delays[i] for i in i_PS]
        self.ps_cluster = {s: x for s, x in zip(good_delay_stations,
                                                good_delays)}

        bad_delay = [s for s in delay_stations if s not in good_delay_stations]
        good_picks = [p for p in picks if p.station not in bad_delay]
        return good_picks

    def _get_origtime_delays(self, vel_mod, target_depth, max_dist=1.,
                             num_vals=20):
        """
        :param max_dist: maximum distance (degrees) at which to calculate
        :param num_vals: number of distances (from 0 to max_distance) at
            which to calculate delats
        :returns: dict with keys 'op' and 'ps' corresponding to sorted lists
            of delay times from origin-to-p and p-to-s
        """
        try:
            model = TauPyModel(model=vel_mod)
        except Exception:
            build_taup_model(vel_mod)   # converts from np or tvel to npz
            model = TauModel.from_file(vel_mod)  # reads npz
        delays = {'op': [0], 'ps': [0]}
        # ps_delays, p_delays = [0], [0]
        step = max_dist / num_vals
        for distance in np.arange(step, max_dist + .001, step):
            arrivals = model.get_travel_times(target_depth, distance,
                                              phase_list=["P", "S"])
            p_delay = [x.time for x in arrivals if x.phase == "P"]
            s_delay = [x.time for x in arrivals if x.phase == "S"]
            ps_delay = s_delay - p_delay
            if len(ps_delay) == 1 and len(p_delay) == 1:
                delays['ps'].append(ps_delay[0])
                delays['op'].append(p_delay[0])
            else:
                log('Could not calculate delays for dist={:.2g} degrees'
                    .format(distance), 'warning')
        return delays

    def _ps_to_otime(self, p_time, s_time):
        """
        Calculate origin time given P and S arrival times

        :param p_time: p arrival time
        :param s_time: s arrival time
        """
        return self.calc_origin_time(p_time, s_time, self.vp_over_vs,
                                     self.delays)

    @staticmethod
    def calc_origin_time(p_time, s_time, vpvs=1.65, delays_model=None):
        """
        Calculate origin time given P and S arrival times

        :param p_time: P arrival time (UTCDateTime)
        :param s_time: S arrival time (UTCDateTime)
        :param delays_model: dict with keys='op' and 'ps', each with an
            equal-length sorted list of origin-to-P and P-to-S delays
        :param vpvs: Vp/Vs ratio, to be used if delays_model == None, using the
            equation  o_time = p_time - ps_delay / (vp/vs - 1)
        :returns: origin_time
        :rtype: UTCDateTime
        """
        origin_time = p_time
        if delays_model is not None:
            origin_time -= np.interp(s_time-p_time,
                                     delays_model['ps'],
                                     delays_model['op'])
        else:
            origin_time -= (s_time - p_time) / (vpvs - 1)
        return origin_time

    def _cluster_clean_otimes(self, times=None):
        """
        Return indices of origin times fitting in largest cluster
        """
        indices = cluster_clean_indices(self.min_clust,
                                        self.cluster_window_otime,  times)
        return [times[i] for i in indices]

    def _cluster_clean_picks(self, phase=None, picks=None):
        """
        Return picks fitting in the largest cluster

        :param phase: "P" or "S"
        :param picks: list of obspy Picks
        :returns: picks that fit the cluster criteria
        """
        assert phase in 'PS', f'phase is "{phase}", should be "P" or "S"'
        if phase == 'P':
            window_sec = self.cluster_window_P
        else:
            window_sec = self.cluster_window_S
        indices = cluster_clean_indices(self.min_clust, window_sec,
                                        [x.time for x in picks])
        return [picks[i] for i in indices]


def cluster_clean_indices(min_clust=None, window_sec=None, times=None):
    """
    Return indices of times fitting in the largest cluster

    :param min_clust: minimum cluster size
    :param window_sec: cluster window length
    :param times: list of UTCDateTime

    :returns: indices of times that fit the cluster criteria
    """
    if len(times) == 0:
        return []
    elif len(times) == 1:
        return [0]
    else:
        cluster_data = fclusterdata(np.array([[x.timestamp for x in times]]).T,
                                    t=window_sec, criterion='distance')
        cluster_groups = dict()
        for cnum in cluster_data:
            cluster_groups[cnum] = cluster_groups.get(cnum, 0) + 1
        max_length = max(cluster_groups.values())
        max_keys = [k for k, v in cluster_groups.items() if v == max_length]

        if len(max_keys) > 1 and max_length >= min_clust:
            log('{:d} clusters of max_length ({:d}), returning first one'
                .format(len(max_keys), max_length), 'warning')
            log(np.nonzero(cluster_data == max_keys[0])[0], 'debug')
        indices = np.nonzero(cluster_data == max_keys[0])[0]
        return indices


def clean_distri(values, n_std=3, mode='median', min_vals=3):
    """
    Remove extreme values in a vector

    :param values: list of values
    :param n_std: erase values for which x is more than n_std standard
        deviations from mode(x)
    :param mode: 'mean' or 'median' [ default = 'median' ]
    :param min_vals:  do not clean if vector is shorter than this[default=3]
    :returns: list of values after cleaning, indices of values retained

    >>> clean_distri([0, 3, 4, 6, 3, 5, 25])
    (array([0, 3, 4, 6, 3, 5]), array([0, 1, 2, 3, 4, 5]))
    """
    if len(values) < min_vals:
        return values, list(range(len(values)))

    x = np.array(values)
    if mode == 'mean':
        dev = x - np.mean(x)
    else:
        dev = x - np.median(x)

    # Calculate std of values
    lim = round(0.68*len(x))  # assume gaussian around median
    if lim == 0:
        lim = 1
    std = np.sort(np.abs(dev))[lim]

    in_rm = np.nonzero(np.abs(dev) <= n_std * std)[0]

    return x[in_rm], in_rm


def _pick_stations(picks):
    # log(picks, 'debug')
    return list(set([x.station for x in picks]))


if __name__ == "__main__":
    import doctest
    doctest.testmod()
