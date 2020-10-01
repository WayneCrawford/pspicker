# Generated with SMOP  0.41-beta
# from smop.libsmop import *

# Standard libraries
import os.path
import warnings
import glob

# Publicly available libraries
import numpy as np
from scipy import stats
from obspy.core import UTCDateTime
from obspy.core.event import Catalog as obspy_Catalog
from obspy.core.inventory.response import PolesZerosResponseStage
from obspy.core.event.base import QuantityError as obspy_QuantityError
from obspy.core.event.magnitude import Amplitude as obspy_Amplitude
from obspy.core.event.origin import Pick as obspy_Pick
from obspy.core.event import Event as obspy_Event
from obspy.core.stream import Stream
from obspy.signal.invsim import simulate_seismometer
from obspy.io.nordic.core import read_nordic
from obspy.core import read as obspy_read

# module libraries
from .parameters import (PickerParameters, PickerRunParameters,
                         PickerLoopParameters)
from .kurtosis import Kurtosis
from .ps_picker_plotter import PSPickerPlotter
from .trace_utils import (fast_polar_analysis, mean_trace, pk2pk, same_inc,
                          select_traces, smooth_filter, snr_function)
from .utils import (clean_distri, cluster_cleaning,
                    get_response, picks_ps_times, picks_matched_stations)


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

#    @property
#    def stations(self):
#        return self.param.station_parameters.stations
#
#    @property
#    def n_stations(self):
#        return self.param.station_parameters.n_stations

    def __str__(self):
        """
        """
        str = "PSPicker\n"
        str += f"    wav_base: {self.wav_base_path}\n"
        str += f"    output db directory: {self.rea_path_out}\n"
        str += f"    parameters: {self.param}\n"
        return str

    def run_many(self, rea_base_path_in, start_yearmonth, end_yearmonth,
                 show_plots=False):
        """
        Loops over events in a date range

        :param rea_base_path_in: absolute basepath to the input datbase
            files (just before the YEAR/MONTH subdirectories)
        :param start_yearmonth: "YYYYMM" of first month to process
        :param end_yearmonth: "YYYYMM" of last month to process
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
                print(f'{year:04d}/{month:02d}')
                s_files = glob.glob(os.path.join(
                    rea_base_path_in, f'{year:04d}', f'{month:02d}', '*.S*'))
                for s_file in s_files:
                    try:
                        self.run_one(s_file, show_plots=show_plots)
                    except Exception as err:
                        warnings.warn('Pick_Function failed for {s_file}')
                        self._write_debug_file(debug_fname, err, s_file)

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

    def run_one(self, rea_name, show_plots=True, debug=False):
        """
        Picks P and S arrivals on one waveform, using the Kurtosis

        Information in the database file will be appended with the picks.
        :param rea_name: database file to read (full path)
        :param show_plots: show plots (useful for training)
        """
        self._setup(rea_name, show_plots)
        self.run.plotter.global_window_timebounds(self.run.global_first_time,
                                                  self.run.global_last_time,
                                                  sorted(self.run.stations))
        self.run.plotter.pick_window_start(self.run.global_first_time,
                                           self.run.global_last_time,
                                           sorted(self.run.stations))

        # Now pick on the traces
        self.debug = debug
        amplitudes = []
        picks = []
        iter = 0
        for station_name, chan_map in self.run.channel_maps.items():
            # Reject stations not listed in parameter file
            if station_name not in self.param.stations:
                continue
            self.loop = PickerLoopParameters(
                station=station_name,
                station_params=self.param.station_parameters[station_name],
                channel_map=chan_map)
            self.loop.add_component_traces(self.run.stream)
            self.run.plotter.station_window_start(self.loop.datP[0], iter)

            onset_P, onset_S = None, None  # samples after self.loop.t_begin

            # SNR analysis
            fe = self.loop.station_params.f_energy
            datS_filtered = self.loop.datS.copy().filter(
                'bandpass', corners=3, freqmin=fe[0], freqmax=fe[1])
            snr, energy = self._calc_snr(datS_filtered)
            if self._snr_trustworthy(snr):
                extrema, onset_P, onset_S = self._run_Kurtosis(energy)

                # Verify phases using Polarity analysis
                if self.loop.station_parms.use_polarity\
                        and (len(datS_filtered) == 3):
                    onset_P, onset_S = self._polarity_analysis(
                        extrema, datS_filtered)

            self.run.plotter.station_window_Ptrace(self, iter)
            self.run.plotter.plot_onsets(self, iter, onset_P, onset_S)

            new_picks = self._make_picks(onset_P, onset_S, snr)
            picks.extend(new_picks)
            amplitudes.append(self._calc_amplitude(new_picks, snr))
            iter += 1

        # Jacknife ###########
        picks = self._remove_unassociated(picks)

        # Only keep biggest amplitude between P and S
        # Replaced by amplitudes.extend above?????
        # G = cells2amp(P_picked_cell, S_picked_cell)

        self.run.plotter.pick_window_add_picks(self, picks)
        self._save_event(picks, amplitudes)

    def _setup(self, rea_name, show_plots):
        """
        Setup starting parameters and objects

        :rea_name: database file to read (full path)
        :show_plots: whether to plot or not
        """
        plotter = PSPickerPlotter(show_plots=show_plots)
        full_wavefile = self._setuprun_get_wavefile_name(rea_name)
        stream = obspy_read(full_wavefile, 'MSEED')
        # get rid of bad last sample in some streams, and detrend
        for tr in stream:
            tr.data = tr.data[:-10]
            tr.detrend(type='demean')
        channel_maps = select_traces(stream, self.param.channel_mapping_rules)
        overall_distri, channel_maps, plotter =\
            self._setuprun_remove_problem_stations(stream, channel_maps,
                                                   plotter)
        ft, lt, overall_distri = self._setuprun_define_analysis_window(
            stream, overall_distri)
        self.run = PickerRunParameters(rea_name=rea_name,
                                       wavefile=full_wavefile,
                                       stream=stream,
                                       channel_maps=channel_maps,
                                       overall_distri=overall_distri,
                                       global_first_time=ft,
                                       global_last_time=lt,
                                       plotter=plotter)

    def _setuprun_get_wavefile_name(self, rea_name):
        """
        Get the full WAV filename

        :param rea_name: full path to database file
        """
        # Pick_Function.m:103
        # full_name = os.path.join(self.rea_path, self.rea_name)
        if self.verbose:
            print(f'database filename = {rea_name}')
        cat, wav_names = read_nordic(rea_name, return_wavnames=True)
        assert len(wav_names) == 1, 'More than one wav_name in database file'
        parts = wav_names[0][0].split('-')
        full_wav_name = os.path.join(self.wav_base_path, parts[0], parts[1],
                                     wav_names[0][0])
        return full_wav_name

    def _setuprun_remove_problem_stations(self, stream, channel_maps, plotter,
                                          n_smooth=15):
        """
        Remove flat-lined stations

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
        plotter.global_window_setup()
        i_trace = 0
        rm_stations = []
        # for iter in arange(1,length(TRZ)).reshape(-1):
        for station, channel_map in sorted(channel_maps.items(),
                                           key=lambda x: x[0]):
            trace = stream.select(id=channel_map.Z)[0]
            sr = trace.stats.sampling_rate
            assert station == trace.stats.station,\
                'trace station ({}) != iteration station ({})'.format(
                    trace.stats.station, station)
            if np.all(np.diff(trace.data) == 0):
                print('Station {} is flat-lined, ignoring'.format(station))
                rm_stations.append(station)
                continue
            kurto_cum, _, _ =\
                Kurtosis.trace2kurto(trace, p.gw_frequency_band,
                                     p.gw_sliding_length, n_smooth)
            if len(kurto_cum.data) == 0:
                warnings.warn('kurto_cum is empty, ignoring station "{}"'
                              .format(station))
                # REMOVE THIS STATION FROM KURTOSIS CALCULATION
                rm_stations.extend(station)
                continue

            mean_kurto = mean_trace(kurto_cum)
            kurto_ext, ind_ext, _ =\
                Kurtosis.follow_extrem(mean_kurto, 'mini', p.gw_n_extrema,
                                       [p.gw_extrema_samples], 'no-normalize',
                                       'no-sense')
            ext_seconds = [i/sr for i in ind_ext]
            if self.debug:
                debug_stream = Stream([trace, mean_kurto, kurto_ext[0]])
                debug_stream[1].stats.channel = 'KUR'
                debug_stream[2].stats.channel = 'EXT'
                print(station, ind_ext, ext_seconds)
                debug_stream.plot(equal_scale=False)
            overall_distri.extend(ext_seconds)
            plotter.global_window_onetrace(trace, i_trace, ext_seconds)
            i_trace += 1

        # REMOVE PROBLEM STATIONS (if necessary)
        channel_maps = {s: v for s, v in channel_maps.items()
                        if s not in rm_stations}
        return overall_distri, channel_maps, plotter

    def _setuprun_define_analysis_window(self, stream, overall_distri):
        """
        Define size of new analysis window

        :param stream: stream containing all data traces
        :param overall_distri: array of extrema on all stations (seconds from t_begin)
        :returns: first_time, last_time, overall_distri
        """
        # Pick_Function.m:204
        n_samples = len(stream[0].data)
        sr = stream[0].stats.sampling_rate
        start_time = stream[0].stats.starttime

        # last_sample is n_samples * pick_window_end in percentage
        last_sample = n_samples * self.param.gw_end_cutoff
        # Cut down picks to those within global bounds
        overall_distri = [s for s in overall_distri if s <= last_sample/sr]

        min_global = center_distri(overall_distri,
                                   self.param.gw_distri_secs)

        T_left = self.param.gw_offsets[0]
        T_right = self.param.gw_offsets[1]
        first_offset = min_global + T_left
        last_offset = min_global + T_right
        if first_offset < 0:
            first_offset = 0
        if last_offset >= n_samples/sr:
            last_offset = n_samples/sr

        return (start_time + first_offset,
                start_time + last_offset,
                overall_distri)

    def _calc_snr(self, datS_filtered):
        """
        Calculate the signal-to-noise relation using the S-wave channel(s)
        """
        # energy = sqrt(sum(filt**2, 2))
        # snr = snr_function(energy, rsample, param.SNR_wind(1),
        #                    param.SNR_wind(2))
        # filt=filterbutter(3, station_param.f_energy(1),
        #                   station_param.f_energy(2), rsample, datS)
        # energy = sqrt(sum(filt**2, 2))
        squared = datS_filtered.copy()
        for tr in squared:
            tr.data = np.power(tr.data, 2)
        stacked = squared.stack()
        assert len(stacked) == 1, 'stacked data has more than one trace!'
        stacked[0].data = np.power(stacked[0].data, 0.5)
        energy = stacked[0]
        snr = snr_function(energy,
                           self.param.SNR_noise_window,
                           self.param.SNR_signal_window)
        return snr, energy

    def _snr_trustworthy(self, snr, n_smooth=100):
        """
        Check if the signal can be trusted or not.

        Considers the signal trustworthy if the SNR crosses the specified
        threshold at least SNR_threshold_crossings times

        :param snr: signal-to-noise ratio trace
        :n_smooth: length of moving average filter to apply before analysis
        """
        # clear('C','ll','rate')
        # stations{iter}
        # Pick_Function.m:380
        snr_smooth = smooth_filter(snr, n_smooth)
        snr_smooth = same_inc(snr_smooth, self.run.global_first_time,
                              self.run.global_last_time)
        SNR_threshold = self.get_SNR_threshold(snr_smooth)
        sign_change = np.diff(np.sign(snr_smooth - SNR_threshold))
        SNR_crossings = len(sign_change[sign_change == 2])
        return SNR_crossings >= self.param.SNR_min_threshold_crossings

    def get_SNR_threshold(self, snr):
        """
        Calculate the signal-to-noise threshold value

        Can be set as an absolute value or as a fraction of the maximum SNR,
        using the parameter SNR_threshold_parameter
        :param snr: trace of smoothed signal-to-noise ratio
        """
        tp = self.param.SNR_threshold_parameter
        if (tp > 0 and tp < 1):
            threshold = tp * snr.data.max()
        else:
            assert tp < 0, f'Illegal SNR_threshold_parameter value: {tp:g}'
            threshold = -tp
        # Minimum possible SNR threshold is the quality = '3' SNR
        if threshold < self.param.SNR_quality_thresholds[0]:
            threshold = self.param.SNR_quality_thresholds[0]
        return threshold

    def _run_Kurtosis(self, snr, energy):
        """
        calculate extrema and estimate onset_P, onset_S using the Kurtosis
        """
        first_time, last_time = self._refine_pick_window(energy)
        all_mean_M, _ = Kurtosis.trace2FWkurto(
            self.loop.datP, self.loop.station_params.Kurto_F,
            self.loop.station_params.Kurto_W, 1, first_time, last_time)
        # Pick_Function.m:430
        kurto_modif, ind_ext, ext_value = Kurtosis.follow_extrem(
            all_mean_M, 'mini', self.loop.station_params.n_follow,
            self.loop.station_params.Kurto_S,
            'no-normalize', 'no-sense')
        extrema = [{'i': i, 'snr': snr[ind_ext]} for i in ind_ext]
        self.run.plotter.station_window_add_snr_nrg(
                self, snr, energy, all_mean_M, kurto_modif)
        if len(extrema) == 0:
            onset_P, onset_S = None, None
        else:
            self.run.plotter.station_window_add_extrema(self, extrema)
            extrema, onset_P, onset_S = self._calc_follows(ind_ext)
        return extrema, onset_P, onset_S

    def _refine_pick_window(self, energy):
        """
        Refine pick window if requested
        """
        # choose the first datP trace as the time and sr reference
        tr = self.loop.datP[0]
        sr = tr.stats.sampling_rate
        tr_start = tr.stats.starttime

        if self.loop.station_params.energy_window == 0:
            return (self.run.global_first_time,
                    self.run.global_last_time)

        energy_smooth = smooth_filter(energy, 50)
        energy_smooth = same_inc(energy_smooth,
                                 self.run.global_first_time,
                                 self.run.global_last_time)
        __, ind_max = np.max(energy_smooth)
        last_sample = ind_max.copy()
        maxKurtoWind = np.max(self.loop.station_params.Kurto_W)
        max_precursor = np.floor(
            sr * (self.loop.station_params.energy_window
                  + maxKurtoWind))
        first_sample = ind_max - max_precursor
        if first_sample < 0:
            first_sample = 0
            last_sample = max_precursor
        # Pick_Function.m:417
        return tr_start + (first_sample/sr, last_sample/sr)

    def _polarity_analysis(self, extrema, datS_filtered):
        """
        Set onset_P and onset_S from extrema based on signal polarity
        """
        # Pick_Function.m:543
        rectP, aziP, dipP = fast_polar_analysis([e['i'] for e in extrema],
                                                2, datS_filtered[0],
                                                datS_filtered[1],
                                                datS_filtered[2])
        # dip_rect = np.multiply(np.sin(np.abs(np.radians(dipP))), rectP)
        dipp = np.sin(np.abs(np.radians(dipP)))
        smooth_dipp = np.lfilter(np.ones(100) / 100, 1, dipp)
        smooth_rectilinP = np.lfilter(np.ones(100) / 100, 1, rectP)
        Drb = np.sign(1.3 * smooth_dipp - smooth_rectilinP)
        DR = np.lfilter(np.ones(200) / 200, 1, np.multiply(rectP, Drb))
        DR[rectP == 0] = 0
        if len(extrema) == 0:
            return None, None
        else:
            onset_P, onset_S = None, None
            if len(extrema) == 1:
                mat_DR = DR(np.arange(np.floor(extrema[0]['i']) - 200,
                                      np.floor(extrema[0]['i']) + 200))
                tmp = np.max(np.abs(mat_DR)) * np.sign(np.mean(mat_DR))
                if tmp >= self.param.dip_rect_thresholds['P']:
                    onset_P = extrema[0]['i']
                elif tmp <= self.param.dip_rect_thresholds['S']:
                    onset_S = extrema[0]['i']
            elif len(extrema) == 2:
                mat_DR1 = DR(np.arange(np.floor(extrema[0]['i']) - 200,
                                       np.floor(extrema[0]['i'] + 200)))
                mat_DR2 = DR(np.arange(np.floor(extrema[1]['i']) - 200,
                                       np.floor(extrema[1]['i'] + 200)))
                if np.mean(mat_DR1) > np.mean(mat_DR2):
                    onset_P = extrema[0]['i']
                    onset_S = extrema[1]['i']
        return onset_P, onset_S

    def _calc_follows(self, extrema, ind_ext):
        """
        returns offsets depending on number of extrema to follow_extrem

        :returns extrema (cleaned), onset_P, onset_S
        :rtype: dict, int, int
        """
        # Pick_Function.m:507
        # eliminate extrema whose snr is less than SNR_thresh
        extrema = [x for x in extrema if x['snr'] > self.param.SNR_thresh]
        assert self.loop.sation_params.nfollow == 1 \
            or self.loop.sation_params.nfollow == 2, \
            'n_follow is not 1 or 2'
        if self.loop.station_params.n_follow == 1:
            if len(extrema) == 0:
                return [], 0, 0
            else:
                return extrema, ind_ext.copy(), 0
        elif self.loop.station_params.n_follow == 2:
            extrema.sort(key='i')
            if len(extrema) == 0:
                return [], 0, 0
            elif len(extrema) == 1:
                return extrema, extrema[0]['i'], 0
            else:
                return extrema, extrema[0]['i'], extrema[1]['i']

    def _make_picks(self, onset_P, onset_S, snr):
        """
        Put onset_P and onset_S into list of obspy Picks
        """
        picks = []
        base_waveid = self.loop.datP[0].get_id()[:-2]
        if onset_P is not None:
            cmp = self.loop.station_params.putPick_P_Comp
            picks.extend(obspy_Pick(
                time=self.loop.index_to_time(onset_P),
                time_errors=self._SNR_to_time_error(snr.data[onset_P]),
                waveform_id=base_waveid + cmp,
                phase_hint=self.loop.station_params.putPick_P_Phase))
        if onset_S is not None:
            cmp = self.loop.station_params.putPick_S_Comp
            picks.extend(obspy_Pick(
                time=self.loop.index_to_time(onset_S),
                time_errors=2*self._SNR_to_time_error(snr.data[onset_S]),
                waveform_id=base_waveid + cmp,
                phase_hint=self.loop.station_params.putPick_S_Phase))
        return picks

    def _calc_amplitude(self, picks, snr):
        """
        Calculate maximum Woods-Anderson amplitude

        :param picks: obspy picks for this station
        :param snr: trace containing signal-to-noise ratio
        :returns: obspy.core.event.magnitude.Amplitude
        """
        pick_P = [x for x in picks if x.phase_hint[0] == 'P']
        pick_S = [x for x in picks if x.phase_hint[0] == 'S']

        # Pick_Function.m:621
        # Look for response file in local directory, then data/CAL directory
        # Should take advantage of obspy to read all standard formats
        filename = self.loop.station_params.response_file
        if not os.path.isfile(filename):
            cal_dir = os.path.join(self.param.Main_path, 'CAL/')
            filename = os.path.join(cal_dir, filename)
        paz = get_response(filename, self.param.response_file_type)
        paz_wa = _wood_anderson_paz()
        # Set gain to 1 so that output will be m/s, not volts
        # paz_wa.stage_gain = 1
        # want output in m?  should I set wa gain and sensitivity to 1?
        if self.loop.station_params.use_polarity:
            wood = self.loop.dat_noH.copy()
        else:
            wood = self.loop.datP.copy()
        for tr in wood:
            tr.simulate(paz_remove=paz, paz_simulate=paz_wa, water_level=60.0)
            #tr.data = simulate_seismometer(
            #    tr.detrend().data, tr.stats.sampling_rate, paz, paz_wa)
        if len(pick_S) > 0:
            pick = pick_S[0]
            Amp = pk2pk(wood, pick.time, before_pick=20, after_pick=10)
        if len(pick_P) > 0:
            pick = pick_P[0]
            Amp = pk2pk(wood, pick.time, before_pick=5, after_pick=30)
        else:
            return None
        amplitude = obspy_Amplitude(
                generic_amplitude=Amp['amplitude'],
                type='AML',
                # category='other', # ["point", "mean", "duration",
                #                      "period", "integral", "other"]
                unit='m/s',  # obspy.core.event.header.AmplitudeUnit,
                period=Amp['period'],
                snr=Amp['snr'],
                pick_id=pick.id,
                waveform_id=pick.waveform_id)
        return amplitude

    def _remove_unassociated(self, picks):
        """
        Very basic (clustering-based) pick selection by association

        DOESN'T add/subsitute in alternative picks for those thrown out
        DOESN'T use a velocity model, max depth and min distance to estimate
            the possible range of pick times
        """
        # Remove lines with arrival times = start_time
        # P_picked_cell=rm_cell_line(P_picked_cell[2] == 0,P_picked_cell)
        # S_picked_cell=rm_cell_line(S_picked_cell[2] == 0,S_picked_cell)

        picks = self._remove_unclustered(picks)
        if self.param.assoc_min_picks <= self.run.n_stations:
            picks = self._remove_badly_distributed(picks)
            picks = self._remove_bad_delays(picks)
        return picks

    def _remove_unclustered(self, picks):
        """
        Remove picks and reassign phases based on clustering

        Very basic association
        :param picks: list of Pick objects
        :returns: modified picks
        """
        # Pick_Function.m:746
        picksP = cluster_cleaning(self.param.assoc_cluster_windows['P'],
                                  [p for p in picks if p.phase_hint[0] == 'P'])
        picksS = cluster_cleaning(self.param.assoc_cluster_windows['S'],
                                  [p for p in picks if p.phase_hint[0] == 'S'])
        return picksP + picksS

    def _remove_badly_distributed(self, picks):
        """
        Remove picks that are well outside of pick distribution

        :param picks: input Picks
        """
        _, iP = clean_distri([p.time for p in picks if p.phase_hint[0] == 'P'],
                             self.param.assoc_min_std, 'median',
                             self.param.assoc_min_picks)
        _, iS = clean_distri([p.time for p in picks if p.phase_hint[0] == 'S'],
                             self.param.assoc_min_std, 'median',
                             self.param.assoc_min_picks)
        return picks[iP] + picks[iS]

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
        _, i_clean = clean_distri(
            delays,
            self.param.assoc_min_std_PtoS,
            'median',
            self.param.assoc_min_picks)
        # Include good delay stations AND non-delay stations
        good_delay = [delay_stations[i] for i in i_clean]
        bad_delay = [s for s in delay_stations if s not in good_delay]
        good_picks = [p for p in picks if p.station not in bad_delay]
        return good_picks

    def _save_event(self, picks, amplitudes):
        """
        Save event to NORDIC file
        """
        # Replaces a large section from Pick_Function.m 840-887
        first_time = UTCDateTime(np.min([p.time.matplotlib_date for p in picks]))
        event = obspy_Event(
            event_type='earthquake',
            picks=picks,
            origin=obspy_Origin(time=estimate_origin_time(picks))
            amplitudes=amplitudes)
        cat = obspy_Catalog(events=[event])
        # How can I change uncertainties to "0", "1", "2", "3"?
        # By creating an associated arrival and setting it's time_weight
        # to the appropriate number (which makes no sense because
        # a weight of zero should have no importance!)
        output_dbfile = os.path.join(self.rea_path_out,
                                     os.path.basename(self.run.rea_name))
        cat.write(output_dbfile,
                  format='NORDIC',
                  evtype='L',
                  wavefiles=[self.run.wavefile],
                  high_accuracy=True)

    def _SNR_to_time_error(self, snr):
        """
        Return approximate pick time errors corresponding to SNR

        :param snr: pick signal-to-noise ratio

        :returns: time_errors
        :rtype: obspy QuantityError
        """
        sr = self.run.datP[0].stats.sampling_rate
        if snr > self.param.SNR_quality_thresholds[3]:
            uncertainty = 2. / sr
        elif snr >= self.param.SNR_quality_thresholds[2]:
            uncertainty = 8. / sr
        elif snr >= self.param.SNR_quality_thresholds[1]:
            uncertainty = 32. / sr
        elif snr >= self.param.SNR_quality_thresholds[0]:
            uncertainty = 128. / sr
        else:
            uncertainty = 2000. / sr
        return obspy_QuantityError(uncertainty)


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
    
    Uses P-S delays if possible
    If not, return earliest pick
    
    :param picks: list of obspy Pick objects
    :param vp_over_vs: assumed velocity ratio
    """
    origin_time = _average_ps_o_time(picks)
    if origin_time is not None:
        return origin_time
    else:
        return sorted(picks, key=lambda p: p.time)[0]


def _average_ps_o_time(picks, vp_over_vs):
    """
    Find origin times for each P-S pair
    
    Uses the equation: o_time = p_time - (s_time - p_time)/(vp/vs - 1)
    :param picks: list of obspy Pick objects
    :param vp_over_vs: assumed velocity ratio
    """
    o_times = []
    for ps_delay, p_time in picks_ps_times(picks):
        o_times.append(p_time - ps_delay/(vp_over_vs - 1))
    if len(o_times) == 0:
        return None
    # Throw out values more than 3 std away
    z = np.abs(starts.zscore(o_times))
    return np.mean(o_times[np.nonzero(z<3)[0]])     


def _wood_anderson_paz():
    """
    Return PoleZerosStage for Wood-Anderson seismometer
    """
    return {'gain': 1.0,
           'poles': [(-6.283-4.7124j), (-6.283+4.7124j)],
           'sensitivity': 2080,
           'zeros': [0 + 0j]}
    # return PolesZerosResponseStage(
    #     pz_transfer_function_type='LAPLACE (RADIANS)',
    #     stage_gain=2800,
    #     stage_gain_frequency=10,
    #     zeros=[0],
    #     poles=[-6.2832-4.7124j, -6.2832+4.7124j],
    #     normalization_frequency=10,
    #     stage_sequence_number=1,
    #     input_units='m/s',
    #     output_units='counts')


if __name__ == '__main__':
    pass
