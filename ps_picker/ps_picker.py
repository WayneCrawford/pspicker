# Generated with SMOP  0.41-beta
# from smop.libsmop import *

# Standard libraries
import os.path
import warnings
import glob
# import logging
# import sys

# Publicly available libraries
import numpy as np
from scipy import stats
from obspy.core import UTCDateTime
from obspy.core.event import Catalog as obspy_Catalog
# from obspy.core.inventory.response import PolesZerosResponseStage
from obspy.core.event.origin import Origin as obspy_Origin
from obspy.core.event import Event as obspy_Event
# from obspy.signal.invsim import simulate_seismometer
from obspy.io.nordic.core import read_nordic
from obspy.core import read as obspy_read

# module libraries
from .parameters import (PickerParameters, PickerRunParameters,
                         PickerStationParameters)
from .kurtosis import Kurtosis
from .energy_snr import EnergySNR
from .polarity import Polarity
from .associator import Associator
from .plotter import Plotter
from .local_amplitude import LocalAmplitude
from .utils import (select_traces, smooth_filter, picks_ps_times)
from .logger import log


class PSPicker():
    """
    Pick P and S arrivals on multiple stations using the Kurtosis and
    do a very basic cluster-based association
    """
    def __init__(self, wav_base_path, parm_file,
                 rea_path_out='./Sfile_directory', verbose=True,
                 debug=False):
        """
        :param wav_base_path: absolute basepath to the waveform files
            (just before the YEAR/MONTH subdirectories)
        :param parm_file: path/name of the parameter file
        :param rea_path_out: path to output REA files
        :param verbose: talk a bit
        """
        self.wav_base_path = wav_base_path
        self.parm_file = parm_file
        self.rea_path_out = rea_path_out
        self.rea_name = None
        self.param = PickerParameters.from_yaml_file(parm_file)
        self.verbose = verbose
        self.debug = debug
        self.run = None

    def __str__(self):
        """
        """
        str = "PSPicker\n"
        str += f"    wav_base: {self.wav_base_path}\n"
        str += f"    output db directory: {self.rea_path_out}\n"
        str += f"    parameters: {self.param}\n"
        return str

    def run_many(self, rea_base_path_in, start_yearmonth, end_yearmonth,
                 plot_global=True, plot_stations=False):
        """
        Loops over events in a date range

        :param rea_base_path_in: absolute basepath to the input datbase
            files (just before the YEAR/MONTH subdirectories)
        :param start_yearmonth: "YYYYMM" of first month to process
        :param end_yearmonth: "YYYYMM" of last month to process
        :param plot_global: show global and overall pick plots
        :param plot_stations: show individual station plots
        """
        self.rea_base_path_in = rea_base_path_in
        start_year = float(start_yearmonth[:4])
        start_month = float(start_yearmonth[4:6])
        end_year = float(end_yearmonth[:4])
        end_month = float(end_yearmonth[4:6])
        debug_fname =\
            'ps_picker_debug_{:04d}.{:02d}-{:04d}.{:02d}_run{}.txt'.format(
                start_year, start_month, end_year, end_month,
                UTCDateTime.now().strftime('%Y.%m.%d'))
        for year in range(start_year, end_year + 1):
            if year == start_year:
                first_month = start_month
            else:
                first_month = 1
            if year == end_year:
                last_month = end_month
            else:
                last_month = 12
            for month in range(first_month, last_month + 1):
                log(f'{year:04d}/{month:02d}')
                s_files = glob.glob(os.path.join(
                    rea_base_path_in, f'{year:04d}', f'{month:02d}', '*.S*'))
                for s_file in s_files:
                    try:
                        self.run_one(s_file, plot_global=plot_global,
                                     plot_stations=plot_stations)
                    except Exception as err:
                        warnings.warn('Pick_Function failed for {s_file}')
                        self._write_debug_file(debug_fname, err, s_file)

    def run_one(self, rea_name, plot_global=True, plot_stations=False,
                assoc=None, verbose=None, debug=False):
        """
        Picks P and S arrivals on one waveform, using the Kurtosis

        Information in the database file will be appended with the picks.
        :param rea_name: database file to read (full path)
        :param plot_global: show global and overall pick plots
        :param plot_stations: show individual station plots
        :param assoc: Associator object (useful for multiple runs with same
            Associator)
        """
        # Run basic Kurtosis and assoc to find most likely window for picks
        self.debug = debug
        if verbose is not None:
            self.verbose = verbose
        plotter = Plotter(plot_global, plot_stations)

        # Read in data and select global pick window
        st, fwf = self._read_waveforms(rea_name)
        cmaps, ft, lt = self._choose_global_window(st, plotter)
        self.run = PickerRunParameters(rea_name=rea_name, wavefile=fwf,
                                       stream=st, channel_maps=cmaps,
                                       first_time=ft, last_time=lt)
        plotter.gw.plot_timebounds(ft, lt, self.run.stations)
        plotter.pw.setup(ft, lt, self.run.stations)

        # Pick on individual traces
        candidates, picks = [], []
        for sta, chan_map in self.run.channel_maps.items():
            # Reject stations not listed in parameter file
            if sta not in self.param.stations:
                if self.verbose:
                    log(f'{sta} not in self.param.stations, ignored',
                        'verbose')
                continue
            p, c = self._pick_one_station(sta, chan_map, plotter)
            picks.extend(p)
            candidates.extend(c)

        if self.verbose:
            log(f'{len(picks):d} picks before association',  level='verbose')
        if assoc is None:
            assoc = Associator(self.param.assoc)
        # picks = assoc.remove_nonclustered(picks)
        picks = assoc.find_same_origin_time(picks, candidates)
        if self.verbose:
            log(f'{len(picks):d} picks after association', level='verbose')
        plotter.pw.plot_picks(picks, self.loop.t_begin, assoc.p_cluster,
                              assoc.s_cluster, assoc.o_cluster)
        obspy_pa = [x.to_obspy(self.run.channel_maps,
                               self.param.SNR.quality_thresholds)
                    for x in picks]
        obspy_picks = [x[0] for x in obspy_pa]
        obspy_arrivals = [x[1] for x in obspy_pa if x[1] is not None]
        amplitudes, obspy_picks = self._calc_amplitudes(obspy_picks)
        self._save_event(obspy_picks, amplitudes, obspy_arrivals)

    def _pick_one_station(self, station_name, chan_map, plotter):
        """
        Calculate picks and candidates for one station

        :param station_name: station name
        :param chan_map: ChannelMap object for the station
        :param plotter: Plotter object
        """
        # make shortened reference to often-used station_parameters
        station_params = self.param.station_parameters[station_name]
        self.loop = PickerStationParameters(station=station_name,
                                            station_params=station_params,
                                            channel_map=chan_map,
                                            stream=self.run.stream)
        plotter.sw.setup(self.loop.datP[0])
        c_P, c_S, candidates = None, None, []

        # SNR analysis
        datS_filt = self.loop.datS.copy().filter(
            'bandpass', corners=3,
            freqmin=station_params.energy_frequency_band[0],
            freqmax=station_params.energy_frequency_band[1])
        energy = EnergySNR(datS_filt, self.param.SNR, plot=self.debug)
        if energy.slice(self.run.first_time,
                        self.run.last_time).is_trustworthy():
            if self.verbose:
                log(f"{station_name}: snr trustworthy", level='verbose')
            c_P, c_S, kurt, candidates = self._run_Kurtosis(energy, plotter)
            for c in candidates:
                c.station = station_name

            # Verify phases using Polarity analysis
            DR = None
            if station_params.use_polarity and (len(datS_filt) == 3):
                c_P, c_S, DR, candidates = self._polarity_analysis(
                        c_P, c_S, candidates, datS_filt)
            plotter.sw.plot_data(self.run.first_time, self.run.last_time,
                                 self.param.SNR.quality_thresholds,
                                 self.loop.datP[0],
                                 energy,
                                 kurt.mean_kurtosis,
                                 kurt.kurto_gradients,
                                 DR,
                                 self.param.polarity.DR_threshold_P,
                                 self.param.polarity.DR_threshold_S)
            if len(candidates) > 0:
                plotter.sw.candidates(candidates,
                                      self.param.polarity.DR_threshold_P,
                                      self.param.polarity.DR_threshold_S)
        elif self.verbose:
            log(f"{station_name}: snr untrustworthy, not picking", 'verbose')

        plotter.pw.plot_traces_candidates(self.loop.datP, c_P, c_S,
                                          candidates, station_name)
        plotter.sw.onsets(c_P, c_S, self.loop.data_limits)

        new_picks = self._make_picks(c_P, c_S)
        return new_picks, candidates

    def _read_waveforms(self, rea_name, format='NORDIC'):
        if format == 'NORDIC':
            full_wavefile = self._get_nordic_wavefile_name(rea_name)
        else:
            raise NameError(f'type {type} not implemented')
        stream = obspy_read(full_wavefile, 'MSEED')
        # get rid of bad last sample in some streams, and detrend
        for tr in stream:
            tr.data = tr.data[:-10]
            tr.detrend(type='demean')
        return stream, full_wavefile

    def _get_nordic_wavefile_name(self, rea_name):
        if self.verbose:
            log(f'database filename = {rea_name}', level='verbose')
        cat, wav_names = read_nordic(rea_name, return_wavnames=True)
        assert len(wav_names) == 1, 'More than one wav_name in database file'
        parts = wav_names[0][0].split('-')
        full_wav_name = os.path.join(self.wav_base_path, parts[0], parts[1],
                                     wav_names[0][0])
        return full_wav_name

    def _choose_global_window(self, stream, plotter):
        """
        Choose the global pick window

        :param stream: all the data read in
        :param plotter: Plotter object
        """
        t_begin = min([t.stats.starttime for t in stream])
        t_end = max([t.stats.endtime for t in stream])
        chan_maps = select_traces(stream, self.param.channel_mapping_rules)
        if self.verbose:
            log(self._channel_maps_str(chan_maps), level='verbose')
        distri, chan_maps = self._gw_get_distri(stream, chan_maps, plotter)
        ft, lt, distri = self._gw_set_window(t_begin, t_end, distri)
        if self.verbose:
            log(f'Global window bounds: {ft} to {lt}', level='verbose')
        return chan_maps, ft, lt

    def _channel_maps_str(self, channel_maps):
        s = 'Channel Map:\n'
        first_iter = True
        for key, v in channel_maps.items():
            if first_iter:
                s += f'{"Station":8s}|{v.__str__(format="table_header")}\n'
                first_iter = False
            s += f'{key:8s}|{v.__str__(format="table_row")}\n'
        return s

    def _gw_get_distri(self, stream, channel_maps, plotter, n_smooth=15):
        """
        Get overall pick distribution (and remove flat-lined stations)

        Finding the Global minimum on all the stations
        Filtering parameters + first window + smooth

        :param stream: all traces
        :param channel_maps: mapping of channel names to components
        :param plotter: the plotter object
        :n_smooth: how many samples to smooth over for calculating Kurtosis
        :returns: overall_distribution of extrema, channel_maps, plotter
        """
        # Pick_Function.m:134
        p = self.param
        overall_distri = []
        i_station = 0
        rm_stations = []
        for station, channel_map in channel_maps.items():
            trace = stream.select(id=channel_map.Z)[0]
            assert station == trace.stats.station,\
                'trace station ({}) != iteration station ({})'.format(
                    trace.stats.station, station)
            if np.all(np.diff(trace.data) == 0):
                log(f'Station {station} flat-lined, ignoring', 'warning')
                rm_stations.append(station)
                continue
            k = Kurtosis([p.gw.kurt_frequency_band], [p.gw.kurt_window_length],
                         n_smooth)
            candidates = k.pick_trace(trace, p.gw.n_extrema,
                                      extrem_smooths=[p.gw.
                                                      kurt_extrema_smoothing])
            for x in candidates:
                x.station = station
            overall_distri.extend([x.time for x in candidates])
            plotter.gw.plot_trace(trace, i_station, candidates)
            i_station += 1

        # REMOVE PROBLEM STATIONS (if necessary)
        channel_maps = {s: v for s, v in channel_maps.items()
                        if s not in rm_stations}
        return overall_distri, channel_maps

    def _gw_set_window(self, t_begin, t_end, overall_distri):
        """
        Define size of new analysis window

        :param t_begin: reference starttime for all traces
        :param t_end: end of all traces
        :param overall_distri: array of extrema on all stations (UTCDateTimes)
        :returns: first_time, last_time, overall_distri
        """
        # Pick_Function.m:204
        # max_offset is data length * gw_end_cutoff
        max_time = t_begin + self.param.gw.end_cutoff * (t_end - t_begin)
        # Cut down picks to those within global bounds
        overall_distri = [t for t in overall_distri if t <= max_time]
        min_global = UTCDateTime(center_distri(
            [t.timestamp for t in overall_distri], self.param.gw.distri_secs))

        first_time = max(min_global + self.param.gw.offsets[0], t_begin)
        last_time = min(min_global + self.param.gw.offsets[1], t_end)
        # if first_time < t_begin:
        #     first_time = t_begin
        # if last_time >= t_end:
        #     last_time = t_end
        return first_time, last_time, overall_distri

    def _run_Kurtosis(self, energy, plotter, debug=False):
        """
        calculate extrema and estimate P and S onsets using the Kurtosis

        :param energy: EnergySNR object
        :param plotter: Plotter object
        :returns:
            c_P: P PickCandidate
            c_S: S PickCandidate
            kurtosis: Kurtosis object
            candidates: list of all PickCandidtes
        """
        first_time, last_time = self._refine_pick_window(energy.nrg)
        if len(self.loop.datP) > 1:
            warnings.warn('Only working on first trace in datP')
        k = Kurtosis(self.loop.station_params.kurt_frequency_bands,
                     self.loop.station_params.kurt_window_lengths,
                     1)
        candidates = k.pick_trace(self.loop.datP[0],
                                  self.loop.station_params.n_extrema,
                                  first_time, last_time,
                                  extrem_smooths=self.loop.station_params.
                                  kurt_extrema_smoothings)
        times = energy.snr.times('utcdatetime')
        for c in candidates:
            c.snr = energy.snr.data[times.searchsorted(c.time)]\

        # calculate picks using only candidates with SNR above min threshold
        strong_candidates = [x for x in candidates
                             if x.snr > min(self.param.SNR.quality_thresholds)]
        c_P, c_S = self._calc_follows(strong_candidates)
        return c_P, c_S, k, candidates

    def _refine_pick_window(self, energy):
        """
        Refine pick window to account for analysis windows?

        Seems like a convoluted way to reduce the time window by
        the longest kurtosis window length + the station energy window?
        I mean, why are we even smoothing the energy window?
        :param energy: energy trace
        """
        if self.loop.station_params.energy_window == 0:
            return (self.run.first_time,
                    self.run.last_time)

        energy_smooth = smooth_filter(energy, 50)
        energy_smooth.trim(self.run.first_time, self.run.last_time)
        starttime = energy_smooth.stats.starttime
        sr = energy_smooth.stats.sampling_rate
        # print(energy_smooth, type(energy_smooth), energy_smooth.data)
        ind_max = np.nanargmax(energy_smooth.data)
        last_sample = ind_max.copy()
        max_kurto_wind = np.max(self.loop.station_params.kurt_window_lengths)
        max_precursor = np.floor(
            sr * (self.loop.station_params.energy_window + max_kurto_wind))
        first_sample = ind_max - max_precursor
        if first_sample < 0:
            first_sample = 0
            last_sample = max_precursor
        # Pick_Function.m:417
        return starttime + first_sample/sr, starttime + last_sample/sr

    def _polarity_analysis(self, c_P, c_S, candidates, datS_filtered):
        """
        Return P and S onsets from candidates based on signal polarity

        :param c_P: P-wave PickCandidate
        :param c_S: S-wave PickCandidate
        :param candidates: list of PickCandidate
        :param datS_filtered: stream of 3 traces used for S picking, filtered
            and with Z last
        :returns:
            c_P: PickCandidate P
            c_S: PickCandidate S
            DR: Dip-rectilinearity trace
            candidates: candidates with DR values added
        """
        if len(candidates) == 0:
            return None, None, None, candidates
        pol = Polarity(datS_filtered, params=self.param.polarity,
                       zcomponents=self.param.channel_mapping_rules.compZ,
                       ncomponents=self.param.channel_mapping_rules.compN,
                       ecomponents=self.param.channel_mapping_rules.compE,
                       verbose=self.verbose)
        DR = pol.calc_dip_rect([c.time for c in candidates])
        if DR is None:
            if self.verbose:
                log("DR not returned, keeping input picks", "verbose")
                return c_P, c_S, DR, candidates

        for c in candidates:
            window = DR.slice(
                c.time, c.time + self.param.polarity.calculate_window).data
            c.DR = np.max(np.abs(window)) * np.sign(np.mean(window))
        c_P, c_S = None, None
        if len(candidates) == 1:
            c = candidates[0]
            if c.DR >= self.param.polarity.DR_threshold_P:
                c_P = c
            elif c.DR <= self.param.polarity.DR_threshold_S:
                c_S = c
        elif len(candidates) == 2:
            if candidates[0].time < candidates[1].time:
                cFirst, cSecond = candidates[0], candidates[1]
            else:
                cFirst, cSecond = candidates[1], candidates[0]
            if cFirst.DR >= self.param.polarity.DR_threshold_P:
                c_P = cFirst
            if cSecond.DR <= self.param.polarity.DR_threshold_S:
                c_S = cSecond
            if c_P is None or c_S is None:
                if cFirst.DR > cSecond.DR:
                    c_P = cFirst
                    c_S = cSecond
        if self.debug:
            log(f'polarity-verified c_P={c_P}, c_S={c_S}', 'debug')
        return c_P, c_S, DR, candidates

    def _calc_follows(self, candidates):
        """
        returns offsets depending on number of extrema to follow

        :param candidates: list of PickCandidate
        :returns candidate_P, candidate_S
        :rtype: PickCandidate, PickCandidate
        """
        # Pick_Function.m:507
        # eliminate extrema whose snr is less than SNR.thresh
        n_extrema = self.loop.station_params.n_extrema
        # assert n_extrema in (1, 2), 'n_extrema is not 1 or 2'
        if len(candidates) == 0:
            return None, None
        if n_extrema == 1:
            return candidates[0], None
        else:
            # log(extrema, level='debug')
            candidates.sort(key=lambda x: x.time)
            if len(candidates) == 1:
                return candidates[0], None
            else:
                return candidates[0], candidates[1]

    def _make_picks(self, candidate_P, candidate_S):
        """
        Put onset_P and onset_S into list of Picks

        :param onset_P: PickCandidate for P
        :param onset_S: PickCandidate for S
        :param snr: signal-to-noise levels Trace
        """
        picks = []
        if candidate_P is not None:
            candidate_P.phase_guess = 'P'
            picks.append(candidate_P)
        if candidate_S is not None:
            candidate_S.phase_guess = 'S'
            picks.append(candidate_S)
        return picks

    def _calc_amplitudes(self, picks):
        amplitudes = []
        stations = list(set([p.waveform_id.station_code for p in picks]))
        for station in stations:
            sta_picks = [p for p in picks
                         if p.waveform_id.station_code == station]
            la = LocalAmplitude(
                self.loop.dat_noH, sta_picks,
                self.param.station_parameters[station].resp_file,
                self.param.response_file_type)
            # log(la, 'debug')
            amp, pick = la.get_iaml(method='wood_calc')
            if amp is not None:
                amplitudes.append(amp)
                picks.append(pick)
        return amplitudes, picks

    @staticmethod
    def save_nordic_event(picks, origin_time, filepath, filename,
                          amplitudes=[], arrivals=[], wavefiles=None,
                          evtype='L', debug=True):
        """
        Save event to NORDIC file

        Arrivals are used to fix the pick weighting: curently direct mapping,
        in version 2.0 will use inverse mapping which correponds to the
        difference between QuakeML and Nordic pick weights
        :param picks: event picks (including amplitude picks)
        :type picks: list
        :param origin_time: event origin time
        :type origin_time: UTCDateTime
        :param amplitudes: event amplitudes
        :type amplitudes: list
        :param arrivals: event arrivals
        :type arrivals: list
        :param wavefiles: list of waveform files corresponding to this event
        :param filepath: path to write to
        :param filename: name of file to write
        :param evtype: event type ('L', 'R' or 'D')
        """
        assert isinstance(picks, list)
        assert isinstance(amplitudes, list)
        assert isinstance(arrivals, list)
        if len(picks) == 0:
            warnings.warn('No picks to save!')
        origin = obspy_Origin(time=origin_time, arrivals=arrivals)
        event = obspy_Event(
            event_type='earthquake',
            picks=picks,
            origins=[origin],
            amplitudes=amplitudes)
        if debug:
            log(event, level='debug')
        cat = obspy_Catalog(events=[event])
        # How to change uncertainties to "0", "1", "2", "3"?
        # By creating an associated arrival and setting it's time_weight
        # to the appropriate number (which makes no sense because
        # a weight of zero should have no importance!)
        output_dbfile = os.path.join(filepath, filename)
        cat.write(output_dbfile, format='NORDIC', evtype=evtype,
                  wavefiles=wavefiles, high_accuracy=True)

    def _save_event(self, picks, amplitudes=[], arrivals=[]):
        """
        Save event to NORDIC file
        """
        # Replaces a large section from Pick_Function.m 840-887
        if len(picks) == 0:
            o_time = self.run.first_time
        else:
            o_time = estimate_origin_time(picks)
        self.save_nordic_event(picks, o_time, self.rea_path_out,
                               os.path.basename(self.run.rea_name),
                               amplitudes=amplitudes,
                               arrivals=arrivals,
                               wavefiles=[self.run.wavefile])

    def _write_debug_file(self, debug_fname, err, s_file):
        """
        Write debugging information after a failed run_one

        :param debug_fname: debug file name
        :param err: returned Exception error
        :param s_file: name of the s_file corresponding to this event
        """
        with open(debug_fname, 'a') as fid:
            fid.write('-'*60 + '\n')
            fid.write(f'{err}\n')
            fid.write('To reproduce the error, type:\n')
            fid.write(f'    picker.run_one("{s_file}"\n\n')


def center_distri(v, win_size, n_steps=1000):
    """
    Return the center of the window containing the most values in an array

    :param v: array of values
    :param win_size: window size
    :param n_steps: number of values between min(v) and max(v) to test
    :returns: center of the distribution
    """
    if len(v) == 0:
        return []

    t = []
    s = []
    v_min, v_max = np.min(v), np.max(v)
    for tv in np.linspace(v_min, v_max, n_steps):
        a = tv - win_size / 2
        b = tv + win_size / 2
        in_range = np.logical_and(v > a, v < b)
        n_in_range = len(np.nonzero(in_range)[0])
        t.append(n_in_range)
        s.append(tv)
    return s[np.argmax(t)]


def estimate_origin_time(picks, vp_over_vs=1.7):
    """
    estimate EQ origin time based on pick times

    Uses P-S delays if possible. If not, return the earliest pick
    :param picks: list of obspy Pick objects
    :param vp_over_vs: assumed velocity ratio
    """
    origin_time = _average_ps_o_time(picks, vp_over_vs)
    if origin_time is not None:
        return origin_time
    else:
        times = list(sorted([p.time for p in picks]))
        # print(f'times={times}, picks={picks}', flush=True)
        return times[0]


def _average_ps_o_time(picks, vp_over_vs):
    """
    Find origin times for each P-S pair

    Uses the equation: o_time = p_time - (s_time - p_time)/(vp/vs - 1)
    :param picks: list of obspy Pick objects
    :param vp_over_vs: assumed velocity ratio
    """
    o_times = []
    ps_delays, p_times = picks_ps_times(picks)
    if ps_delays is None:
        return None
    for ps, p in zip(ps_delays, p_times):
        o_times.append(p - ps/(vp_over_vs - 1))
    # Throw out values more than 3 std away
    zs = np.abs(stats.zscore([x.timestamp for x in o_times]))
    mean_timestamp = np.mean([o.timestamp for o, z in zip(o_times, zs)
                              if z < 3])
    return UTCDateTime(mean_timestamp)


if __name__ == '__main__':
    pass
