# Generated with SMOP  0.41-beta
# from smop.libsmop import *

# Standard libraries
from pathlib import Path
import shutil
import warnings
# import glob
# import warnings
# from logging import info
from datetime import datetime
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
from .pick_candidate import PickCandidate
from .associator import Associator
from .plotter import Plotter
from .local_amplitude import LocalAmplitude
from .utils import (select_traces, smooth_filter, picks_ps_times)
from .logger import setup_log, log
from .timer import Timer

warnings.filterwarnings("ignore",
                        message="Lines of type I have not been implemented "
                                "yet, please submit a development request")


class PSPicker():
    """
    Pick P and S arrivals on multiple stations using the Kurtosis and
    do a very basic cluster-based association
    """
    def __init__(self, parm_file, wav_base_path, database_path_in,
                 database_path_out='./Sfile_directory',
                 database_format='NORDIC'):
        """
        :param parm_file: path/name of the parameter file
        :param wav_base_path: absolute basepath to the waveform files
            (just before the YEAR/MONTH subdirectories)
        :param database_path_in: absolute basepath to the database/catalog
            file(s) (just before the YEAR/MONTH subdirectories)
        :param database_path_out: path to output database files
        :param database_format:
            'NORDIC': assume waveform files and database files are named using
                SEISAN conventions and located in YEAR/MONTH subdirectories
                under wav_base_path and database_path_in, respectively
        """
        self.parm_file = parm_file
        self.wav_base_path = Path(wav_base_path)
        self.database_path_in = Path(database_path_in)
        self.database_path_out = Path(database_path_out)
        if not self.database_path_out.is_dir():
            assert not self.database_path_out.exists()
            self.database_path_out.mkdir()
        # self.database_filename = None
        self.param = PickerParameters.from_yaml_file(parm_file)
        self.plot_debug = False
        self.run = None
        self.assoc = None
        # self.log_level = None

    def __str__(self):
        """
        """
        str = "PSPicker\n"
        str += f"    wav_base: {self.wav_base_path}\n"
        str += f"    output db directory: {self.database_path_out}\n"
        str += f"    parameters: {self.param}\n"
        return str

    def run_many(self, start_date, end_date, plot_global=False,
                 plot_stations=False, ignore_fails=True, log_level='info'):
        """
        Loops over events in a date range

        :param start_date: first data to process
        :param end_date: last data to process
        :type start_date, end_date: str of format "YYYYMMDD" or "YYYYMMDDHHMM"
        :param plot_global: show global and overall pick plots
        :param plot_stations: show individual station plots
        :param ignore_fails: keep going if one run fails
        :param log_level: console log level (choices = 'debug', 'verbose',
            'info', 'warning', 'error', 'critical'), default='info'
        """
        setup_log(log_level)
        # self.log_level = log_level
        start_dt = self._split_date(start_date)
        end_dt = self._split_date(end_date)
        debug_fname = 'ps_picker_debug_{}-{}_run{}.txt'.format(
                start_date, end_date, datetime.today().strftime('%Y%m%d'))
        log('Running from {} to {}'.format(start_dt.strftime("%Y%m%d-%H%M"),
                                           end_dt.strftime("%Y%m%d-%H%M")))

        # Print parameter information
        log(str(self.param), 'verbose')

        def _date_match(y, m, d, ref):
            return y == ref.year and m == ref.month and d == ref.day

        for year in range(start_dt.year, end_dt.year + 1):
            first_month = 1
            last_month = 12
            if year == start_dt.year:
                first_month = start_dt.month
            if year == end_dt.year:
                last_month = end_dt.month
            for month in range(first_month, last_month + 1):
                first_day = 1
                last_day = 31
                if year == start_dt.year and month == start_dt.month:
                    first_day = start_dt.day
                    # first_hour = start_dt.hour
                    # first_minute = start_dt.minute
                if year == end_dt.year and month == end_dt.month:
                    last_day = end_dt.day
                for day in range(first_day, last_day + 1):
                    kwargs = {}
                    if _date_match(year, month, day, start_dt):
                        kwargs['first_hour'] = start_dt.hour
                        kwargs['first_minute'] = start_dt.minute
                    if _date_match(year, month, day, end_dt):
                        kwargs['last_hour'] = end_dt.hour
                        kwargs['last_minute'] = end_dt.minute
                    self._run_one_day(year, month, day, plot_global,
                                      plot_stations, ignore_fails,
                                      debug_fname, **kwargs)

    def run_one(self, database_filename, plot_global=True, plot_stations=False,
                assoc=None, log_level='verbose', plot_debug=None):
        """
        Picks P and S arrivals on one waveform, using the Kurtosis

        Information in the database file will be appended with the picks.
        :param database_filename: database file to read
        :param plot_global: show global and overall pick plots
        :param plot_stations: show individual station plots
        :param assoc: Associator object (used for multiple runs with same
            Associator)
        :param log_level: console log level (choices = 'debug', 'verbose',
            'info', 'warning', 'error', 'critical').  If None, do not setup log
        :param plot_debug: plot some "debugging" plots
        """
        # if not self.log_level:
        #     setup_log(log_level)
        #     self.log_level = log_level
        if log_level is not None:
            setup_log(log_level)
        if self.assoc is None:
            self.assoc = Associator(self.param.assoc)
        if plot_debug is not None:
            self.plot_debug = plot_debug
        log(f'running {database_filename}', 'debug')
        timer = Timer(logger=None)
        timer.start()
        # Run basic Kurtosis and assoc to find most likely window for picks
        # with Timer(text="Setup Plotter: {:0.4f}s"):
        plotter = Plotter(plot_global, plot_stations)

        # Read in data and select global pick window
        # with Timer(text="Read waveforms: {:0.4f}s"):
        st, wavefile = self._read_waveforms(
            self._full_nordic_database_filename(database_filename))
        # print(st.__str__(extended=True))
        if len(st) == 0:
            log('No data found in {wavefile}, referred by {database_filename}',
                'error')
            return
        sta_list = sorted(list(set([tr.stats.station for tr in st])))
        log('Read waveforms from stations {}'.format(', '.join(sta_list)),
            'verbose')
        # with Timer(text="Choose global window: {:0.4f}s"):
        cmaps, ft, lt = self._choose_global_window(st, plotter)
        _check_timelimits(st, ft, lt)
        self.run = PickerRunParameters(
            database_filename=database_filename, wavefile=wavefile,
            stream=st, channel_maps=cmaps, first_time=ft, last_time=lt)
        plotter.gw.plot_pickbounds(ft, lt)
        plotter.pw.setup(ft, lt, self.run.stations)

        # Pick on individual traces
        candidates, picks = [], []
        for sta, chan_map in self.run.channel_maps.items():
            # Reject stations not listed in parameter file
            if sta not in self.param.stations:
                log(f'{sta} not in self.param.stations, ignored', 'warning')
                continue
            # with Timer(text="Pick one station: {:0.4f}s"):
            p, c = self._pick_one_station(sta, chan_map, plotter)
            picks.extend(p)
            candidates.extend(c)

        # with Timer(text="Associate: {:0.4f}s"):
        picks = self.assoc.run(picks, candidates)
        plotter.pw.plot_picks(picks, self.loop.t_begin, self.assoc)
        # with Timer(text="Save picks: {:0.4f}s"):
        # log(f'picks = {picks}', 'debug')
        picks = PickCandidate.remove_duplicates(picks)
        obspy_pa = [x.to_obspy(self.run.channel_maps,
                               self.param.SNR.quality_thresholds)
                    for x in picks]
        obspy_picks = [x[0] for x in obspy_pa]
        obspy_arrivals = [x[1] for x in obspy_pa if x[1] is not None]
        # amplitudes, obspy_picks = self._calc_amplitudes(obspy_picks)
        amplitudes, amp_picks = self._calc_amplitudes(obspy_picks)
        obspy_picks.extend(amp_picks)
        self._save_event(obspy_picks, amplitudes, obspy_arrivals)
        elapsed_time = timer.stop()
        log('    {}: {:2d} Picks and {:2d} Amplitudes on {:2d} stations in '
            '{:0.2f} seconds'.format(database_filename,
                                     len(obspy_picks) - len(amplitudes),
                                     len(amplitudes), len(cmaps),
                                     elapsed_time))

    def _run_one_day(self, year, month, day, plot_global, plot_stations,
                     ignore_fails, debug_fname, first_hour=None,
                     first_minute=None, last_hour=None, last_minute=None):
        """Select and run events for one day"""
        log(f'Running {year}-{month}-{day}, {first_hour=}, '
            f'{first_minute=}, {last_hour=}, {last_minute=}', 'debug')
        db_path_in = self.database_path_in / f'{year:04d}' / f'{month:02d}'
        s_files = list(db_path_in.glob(f'{day:02d}-*.S*'))
        if len(s_files) > 0:
            if first_hour is not None or first_minute is not None:
                s_files = [f for f in s_files if self._nordic_fname_after(
                           f.name, first_hour, first_minute)]
            if last_hour is not None or last_minute is not None:
                s_files = [f for f in s_files if self._nordic_fname_before(
                           f.name, last_hour, last_minute)]
        s_files.sort()
        if len(s_files) > 0:
            log('Running {:d} events on {:04d}-{:02d}-{:02d}'.format(
                len(s_files), year, month, day))
            for s_file in s_files:
                log("   Running {}...".format(s_file), 'verbose')
                t = Timer(logger=None)
                t.start()
                try:
                    self.run_one(s_file, plot_global=plot_global,
                                 plot_stations=plot_stations, log_level=None)
                except Exception as err:
                    log(f'run_one() failed for {s_file}', 'critical')
                    log(err, 'error')
                    if not ignore_fails:
                        raise Exception(err)
                    log('copying original s-file to dest', 'info')
                    inf = s_file
                    outf = self.database_path_out / Path(s_file).name
                    log(f'{inf} to {outf}')
                    shutil.copyfile(inf, outf)
                # elapsed_time = t.stop()

    @staticmethod
    def _nordic_fname_after(f, hour, minute):
        """
        Returns True if the event is on or after the given hour-minute
        :param f: filename (str)
        :param first_hour:
        :param first_minute:
        """
        ev_hour, ev_minute = int(f[3:5]), int(f[5:7])
        if ev_hour > hour:
            return True
        elif ev_hour == hour and ev_minute >= minute:
            return True
        return False

    @staticmethod
    def _nordic_fname_before(f, hour, minute):
        """
        Returns True if the event is on or before the given hour-minute
        :param f: filename (str)
        :param hour:
        :param minutes:
        """
        ev_hour, ev_minute = int(f[3:5]), int(f[5:7])
        if ev_hour < hour:
            return True
        elif ev_hour == hour and ev_minute <= minute:
            return True
        return False

    @staticmethod
    def _split_date(the_date):
        """Splits a date string ("YYYYMMDD" or "YYMMDDHHMM)
        :returns: datetime object year, month, day, hour, minute
        """
        d = the_date
        assert isinstance(d, str), 'date is not a string'
        assert d.isnumeric(), 'date is not numeric'
        assert len(d) == 8 or len(d) == 12,\
            f'"date" ({d}) is {len(d)} numerals, should be 8 or 12'
        if len(d) == 8:
            return datetime(int(d[:4]), int(d[4:6]), int(d[6:]))
        if len(d) == 12:
            return datetime(int(d[:4]), int(d[4:6]), int(d[6:8]), int(d[8:10]),
                            int(d[10:]))

    def _pick_one_station(self, station_name, chan_map, plotter):
        """
        Calculate picks and candidates for one station

        :param station_name: station name
        :param chan_map: ChannelMap object for the station
        :param plotter: Plotter object
        """
        # make shortened reference to often-used station_parameters
        # with Timer(text="  pick_one_station(): setup {:0.4f}s"):
        station_params = self.param.station_parameters[station_name]
        self.loop = PickerStationParameters(station=station_name,
                                            station_params=station_params,
                                            channel_map=chan_map,
                                            stream=self.run.stream)
        plotter.sw.setup(self.loop.datP[0])
        c_P, c_S, candidates = None, None, []

        # SNR analysis
        # with Timer(text="  pick_one_station(): SNR {:0.4f}s"):
        datS_filt = self.loop.datS.copy().filter(
            'bandpass', corners=3,
            freqmin=station_params.SNR_energy.frequency_band[0],
            freqmax=station_params.SNR_energy.frequency_band[1])
        energy = EnergySNR(datS_filt, self.param.SNR, plot=self.plot_debug)
        trust, message = energy.slice(self.run.first_time,
                                      self.run.last_time).is_trustworthy()
        log(f"{station_name}: SNR {message}", 'verbose')
        if trust:
            # with Timer(text="  pick_one_station(): Kurtosis {:0.4f}s"):
            c_P, c_S, kurt, candidates = self._run_Kurtosis(energy, plotter)
            for c in candidates:
                c.station = station_name

            # Verify phases using Polarity analysis
            DR = None
            # with Timer(text="  pick_one_station(): Polarity analysis
            #      {:0.4f}s"):
            if station_params.use_polarity and (len(datS_filt) == 3):
                c_P, c_S, DR, candidates = self._polarity_analysis(
                        c_P, c_S, candidates, datS_filt)
            plotter.sw.plot_data(self.run.first_time, self.run.last_time,
                                 self.param.SNR.quality_thresholds,
                                 self.loop.datP[0],
                                 energy,
                                 kurt,
                                 DR,
                                 self.param.polarity.DR_threshold_P,
                                 self.param.polarity.DR_threshold_S)
            if len(candidates) > 0:
                plotter.sw.candidates(candidates,
                                      self.param.polarity.DR_threshold_P,
                                      self.param.polarity.DR_threshold_S)

        plotter.pw.plot_traces_candidates(self.loop.datP, c_P, c_S,
                                          candidates, station_name, self.assoc)
        plotter.sw.onsets(c_P, c_S, self.loop.data_limits)

        new_picks = self._make_picks(c_P, c_S)
        return new_picks, candidates

    def _read_waveforms(self, database_filename, format='NORDIC'):
        if format == 'NORDIC':
            full_wavefile = self._get_nordic_wavefile_name(database_filename)
        else:
            raise NameError(f'type {type} not implemented')
        stream = obspy_read(full_wavefile, 'MSEED')
        # get rid of bad last sample in some streams, and detrend
        for tr in stream:
            tr.data = tr.data[:-10]
            tr.detrend(type='demean')
        return stream, full_wavefile

    def _get_nordic_wavefile_name(self, database_filename):
        log(f'database filename = {database_filename}', 'verbose')
        cat, wav_names = read_nordic(database_filename, return_wavnames=True)
        assert len(wav_names) == 1, 'More than one wav_name in database file'
        pts = wav_names[0][0].split('-')
        full_wav_name = self.wav_base_path / pts[0] / pts[1] / wav_names[0][0]
        return str(full_wav_name)

    def _full_nordic_database_filename(self, filename):
        """
        Look for a NORDIC database file and return its path

        Looks in local directory, then self.database_path_in, then
        self.database_path_in/YEAR/month

        :param filename: name of the file to look FormatControl
        :type filename: str or pathlib.Path
        """
        # In local directory
        if Path(filename).is_file():
            return filename
        # In database directory
        fullname = self.database_path_in / filename
        if fullname.is_file():
            return fullname
        # In database/year/month directory
        a = filename.split('.')[-1]
        year = a[1:5]
        month = a[5:]
        fullname = self.database_path_in / year / month / filename
        if fullname.is_file():
            return fullname
        else:
            raise NameError(f'database file "{filename}" not found')

    def _choose_global_window(self, stream, plotter):
        """
        Choose the global pick window

        :param stream: all the data read in
        :param plotter: Plotter object
        """
        t_begin = min([t.stats.starttime for t in stream])
        t_end = max([t.stats.endtime for t in stream])
        chan_maps = select_traces(stream, self.param.channel_mapping_rules)
        plotter.gw.setup(t_begin, t_end, [s for s in chan_maps.keys()])
        log(self._channel_maps_str(chan_maps), 'verbose')
        distri, chan_maps = self._gw_get_distri(stream, chan_maps, plotter)
        ft, lt, distri = self._gw_set_window(t_begin, t_end, distri)
        log(f'Global window bounds: {ft} to {lt}', 'verbose')
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
        :n_smooth: samples to smooth kurtosis over
        :returns: overall_distribution of extrema, channel_maps, plotter
        """
        # Pick_Function.m:134
        p = self.param
        overall_distri = []
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
            p.gw.kurtosis.n_smooth = n_smooth
            k = Kurtosis(p.gw.kurtosis)
            candidates = k.pick_trace(trace, p.gw.max_candidates)
            for x in candidates:
                x.station = station
            overall_distri.extend([x.time for x in candidates])
            plotter.gw.plot_trace(trace, station, candidates)

        # REMOVE PROBLEM STATIONS (if necessary)
        channel_maps = {s: v for s, v in channel_maps.items()
                        if s not in rm_stations}
        log(f'all global picks: {overall_distri}', 'debug')
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
        k = Kurtosis(self.loop.station_params.kurtosis)
        candidates = k.pick_trace(self.loop.datP[0],
                                  self.loop.station_params.max_candidates,
                                  first_time, last_time)
        #  Trace.times('utcdatetime') takes 0.3s per call!
        times = energy.snr.times('timestamp')
        imax = energy.snr.stats.npts - 1
        for c in candidates:
            c.snr = energy.snr.data[min(times.searchsorted(c.time.timestamp),
                                        imax)]

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
        if self.loop.station_params.SNR_energy.window == 0:
            return (self.run.first_time,
                    self.run.last_time)

        energy_smooth = smooth_filter(energy, 50)
        energy_smooth.trim(self.run.first_time, self.run.last_time)
        starttime = energy_smooth.stats.starttime
        sr = energy_smooth.stats.sampling_rate
        # print(energy_smooth, type(energy_smooth), energy_smooth.data)
        ind_max = np.nanargmax(energy_smooth.data)
        last_sample = ind_max.copy()
        max_kurto_wind = np.max(self.loop.station_params.kurtosis
                                                        .window_lengths)
        max_precursor = np.floor(
            sr * (self.loop.station_params.SNR_energy.window + max_kurto_wind))
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
                       zcomponents=self.param.channel_mapping_rules.component_orientation_codes.Z,
                       ncomponents=self.param.channel_mapping_rules.component_orientation_codes.N,
                       ecomponents=self.param.channel_mapping_rules.component_orientation_codes.E)
        DR = pol.calc_dip_rect([c.time for c in candidates])
        if DR is None:
            log("DR not returned, keeping input picks", "debug")
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
        log(f'polarity-verified c_P={c_P}, c_S={c_S}', 'verbose')
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
        max_candidates = self.loop.station_params.max_candidates
        # assert max_candidates in (1, 2), 'max_candidates is not 1 or 2'
        if len(candidates) == 0:
            return None, None
        if max_candidates == 1:
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
        amp_picks = []
        stations = list(set([p.waveform_id.station_code for p in picks]))
        for station in stations:
            sta_picks = [p.copy() for p in picks
                         if p.waveform_id.station_code == station]
            # Inefficient way to get dat_noH
            temp = PickerStationParameters(
                station=station,
                station_params=self.param.station_parameters[station],
                channel_map=self.run.channel_maps[station],
                stream=self.run.stream)
            la = LocalAmplitude(
                temp.dat_noH, sta_picks,
                self.param.station_parameters[station].resp_file,
                self.param.response_file_type)
            # log(la, 'debug')
            amp, pick = la.get_iaml(method='wood_calc')
            if amp is not None:
                amplitudes.append(amp)
                amp_picks.append(pick)
        return amplitudes, amp_picks

    @staticmethod
    def save_nordic_event(picks, origin_time, filepath, filename,
                          amplitudes=[], arrivals=[], wavefiles=None,
                          evtype='L', debug=False):
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
        # if len(picks) == 0:
        #     log('No picks to save!', 'warning')
        origin = obspy_Origin(time=origin_time, arrivals=arrivals)
        event = obspy_Event(
            event_type='earthquake',
            picks=picks,
            origins=[origin],
            amplitudes=amplitudes)
        log(event, 'verbose')
        cat = obspy_Catalog(events=[event])
        # How to change uncertainties to "0", "1", "2", "3"?
        # By creating an associated arrival and setting it's time_weight
        # to the appropriate number (which makes no sense because
        # a weight of zero should have no importance!)
        output_dbfile = Path(filepath) / filename
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
        self.save_nordic_event(picks, o_time, self.database_path_out,
                               Path(self.run.database_filename).name,
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
            fid.write(f'    picker.run_one("{s_file}")\n\n')


def _check_timelimits(st, ft, lt):
    """ Check if there are traces outside of the time limits """
    bad_traces = [tr for tr in st
                  if tr.stats.starttime > lt or tr.stats.endtime < ft]
    if len(bad_traces) > 0:
        stream = st.copy()
        stream.traces = bad_traces
        log('some traces are outside of the pick window:', 'error')
        log(stream, 'error')


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


def _average_ps_o_time(picks, vp_vs):
    """
    Find origin times for each P-S pair

    :param picks: list of obspy Pick objects
    :param vp_vs: assumed velocity ratio
    """
    o_times = []
    ps_delays, p_times = picks_ps_times(picks)
    if ps_delays is None:
        return None
    for ps, p in zip(ps_delays, p_times):
        o_times.append(Associator.calc_origin_time(p, p+ps, vpvs=vp_vs))
    # Throw out nans
    o_ts = [o.timestamp for o in o_times]
    # log(o_ts, 'debug')
    if np.any(np.isnan(o_ts)):
        log('there are o_times with timestamp = NaN!', 'error')
        o_ts = [o for o in o_ts if not np.isnan(o)]
    if len(o_ts) == 1:
        mean_timestamp = o_ts[0]
    elif len(o_ts) == 0:
        log('No valid origin times', 'warning')
        return None
    else:
        # Throw out values more than 3 std away
        zs = np.abs(stats.zscore([x for x in o_ts]))
        # log(zs, 'debug')
        if np.any(np.isnan(zs)):
            log("can't use z-values because they have NaNs "
                f"({len(o_ts)} origin_time(s))", 'warning')
            mean_timestamp = np.mean(o_ts)
        else:
            mean_timestamp = np.mean([o for o, z in zip(o_ts, zs) if z < 3])
    # log(f'mean_timestamp = {mean_timestamp}', 'debug')
    return UTCDateTime(mean_timestamp)


if __name__ == '__main__':
    pass
