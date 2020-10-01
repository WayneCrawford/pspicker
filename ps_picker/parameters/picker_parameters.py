# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# readmain.m
import yaml

from .channel_mapping_rules import ChannelMappingRules
from .station_parameters import StationParameters
# from .other_parameters import OtherParameters


class PickerParameters():
    """
    Picker Parameters
    """
    def __init__(self,
                 gw_frequency_band,
                 gw_sliding_length,
                 gw_offsets,
                 SNR_signal_window,
                 SNR_noise_window,
                 SNR_quality_thresholds,
                 SNR_min_threshold_crossings,
                 assoc_cluster_windows,
                 station_parameters,
                 channel_mapping_rules,
                 dip_rect_thresholds={'P':-.4, 'S':-0.4},
                 gw_distri_secs=2,
                 gw_n_extrema=5,
                 gw_extrema_samples=40,
                 gw_end_cutoff=0.9,
                 SNR_threshold_parameter=0.2,
                 assoc_min_picks=1e9,
                 assoc_max_std=3.2,
                 assoc_max_std_PtoS=4,
                 response_file_type=''):
        """
        Initialize Picker Parameters

        :param gw_frequency_band: The frequency band for kurtosis calculation
            during global rewindowing  [low, high]
        :param gw_sliding_length: Length in seconds of sliding window for
            Kurtosis during global rewindowing
        :param gw_distri_secs: Length in seconds of sliding window to use
            to find the densest pick time over all stations
        :param gw_n_extrema: Number of extrema to use for each trace
            when looking at overall extrema distribution
        :param gw_extrema_samples: Number of samples to use in smoothing
            window when calculating extrema
        :param gw_offsets: Offsets in seconds from densest pick
            time for global (all station based) rewindowing [left, right]
        :param gw_end_cutoff: What fraction of data (from start) to
            look at for global Kurtosis window. 1.0 looks everywhere
        :param SNR_noise_window: Length of SNR noise window in seconds
        :param SNR_signal_window: Length of SNR signal window in seconds
        :param SNR_quality_thresholds: SNR thresholds for pick quality
            list with values for quality = '3', '2', '1' and '0', in order
            the minimum value is also used as a baseline for SNR_threshold
        :param SNR_threshold_parameter: Threshold for SNR-based quality
            evaluation.
            if > 0 and less than 1, then SNR_threshold is this times the
                maximum SNR
            if < 0, then SNR_threshold is abs(SNR_threshold_parameter)
        :param SNR_min_threshold_crossings: Minimum crossings of SNR to accept
            a trace
        :param dip_rect_thresholds: DipRect thresholds {'P': value, 'S': value}
        :param assoc_cluster_windows: Windows in seconds for clustering rejection
            {'P': value, 'S':value}
        :param station_parameters: dict with key=station names,
            items=StationParameters object
        :param channel_mapping_rules: ChannelMappingRules object
        :param assoc_min_picks: Minimum number of values (P-picks, S-picks
            or PS delays) needed to perform distribution-based rejection
            (if > n_stations, will not attempt distribution-based rejection
        :param assoc_max_std: maximum number of deviations from standard
            distribution to accept for P picks, and for S picks
        :param cleandistri_PtoS: maximum number of deviations from standard
            distribution to accept for P-S delays
        :param response_file_type: Format of the response file(s).  'GSE' or
            empty.  If empty use Baillardesque PoleZero format.
        """
        sqt = SNR_quality_thresholds
        sqt_txt = 'SNR_quality_thresholds'
        assert len(sqt) == 4, f'Need 4 {sqt_txt}, found {len(sqt):d}'
        assert sqt == sorted(sqt), f"{sqt_txt} not increasing: {sqt}"
        assert len(gw_frequency_band) == 2, "len(gw_frequency_band) != 2"
        assert len(gw_offsets) == 2, "len(gw_offsets) != 2"

        self.gw_frequency_band = gw_frequency_band
        self.gw_sliding_length = gw_sliding_length
        self.gw_offsets = gw_offsets
        self.gw_end_cutoff = gw_end_cutoff
        self.gw_distri_secs = gw_distri_secs
        self.gw_n_extrema = gw_n_extrema
        self.gw_extrema_samples = gw_extrema_samples
        self.SNR_signal_window = SNR_signal_window
        self.SNR_noise_window = SNR_noise_window
        self.SNR_quality_thresholds = SNR_quality_thresholds
        self.SNR_min_threshold_crossings = SNR_min_threshold_crossings
        self.dip_rect_thresholds = dip_rect_thresholds
        self.assoc_cluster_windows = assoc_cluster_windows
        self.station_parameters = station_parameters
        self.channel_mapping_rules = channel_mapping_rules
        self.SNR_threshold_parameter = SNR_threshold_parameter
        self.assoc_min_picks = float(assoc_min_picks)
        self.assoc_max_std = assoc_max_std
        self.assoc_max_std_PtoS = assoc_max_std_PtoS
        self.response_file_type = response_file_type

    # def _make_seisan_paths(self, name):
    #     """
    #     make seisan paths covering start_date to end_date
    #     :param name: name of the subdirectory ('REA' or 'WAV')
    #     """
    #     paths = []
    #     year = self.start_year
    #     month = self.start_month
    #     while year <= self.end_year:
    #         if year == self.end_year and month > self.end_month:
    #             break
    #         if month > 12:
    #             year += 1
    #             month = 1
    #             continue
    #         paths.extend(os.path.join(self.input_path, 'REA', self.base_name,
    #                                   f'{year:04d}', f'{month:02d}', ''))
    #         month += 1
    #     return paths

    @property
    def stations(self):
        return [s for s in self.station_parameters.keys()]

    @property
    def n_stations(self):
        return len(self.station_parameters)

    def __str__(self):
        str = "PickerParamters:\n"
        str += f"    gw_frequency_band = {self.gw_frequency_band}\n"
        str += f"    gw_sliding_length = {self.gw_sliding_length}\n"
        str += f"    gw_distri_secs = {self.gw_distri_secs}\n"
        str += f"    gw_offsets = {self.gw_offsets}\n"
        str += f"    gw_end_cutoff = {self.gw_end_cutoff}\n"
        str += f"    gw_n_extrema = {self.gw_n_extrema}\n"
        str += f"    gw_extrema_samples = {self.gw_extrema_samples}\n"
        str += f"    SNR_signal_window = {self.SNR_signal_window}\n"
        str += f"    SNR_quality_thresholds = "
        str += f"{self.SNR_quality_thresholds}\n"
        str += f"    SNR_min_threshold_crossings = "
        str += f"{self.SNR_min_threshold_crossings}\n"
        str += f"    dip_rect_thresholds = {self.dip_rect_thresholds}\n"
        str += f"    assoc_cluster_windows = {self.assoc_cluster_windows}\n"
        str += f"    stations = {','.join(self.stations)}\n"
        str += f"    SNR_threshold_parameter = "
        str += f"{self.SNR_threshold_parameter}\n"
        str += f"    assoc_min_picks = {self.assoc_min_picks}\n"
        str += f"    assoc_max_std = {self.assoc_max_std}\n"
        str += f"    assoc_max_std_PtoS = {self.assoc_max_std_PtoS}\n"
        str += f"    response_file_type = '{self.response_file_type}'\n"
        str += f"    channel_mapping_rules = {self.channel_mapping_rules}\n"
        return str

    @classmethod
    def from_yaml_file(cls, filename):
        """
        Read in parameters from a YAML file
        """
        with open(filename, 'r') as fic:
            params = yaml.load(fic)

        # Fill in station parameters
        sp = dict()
        kurtosis = params['kurtosis']
        for station, values in params['station_parameters'].items():
            k_parms = values['k_parms']
            sp[station] = StationParameters(
                P_comp=values['P_comp'],
                S_comp=values['S_comp'],
                f_energy=values['f_nrg'],
                frequencies=kurtosis['frequency_bands'][k_parms['freqs']],
                window_lengths=kurtosis['window_lengths'][k_parms['wind']],
                smoothing_sequence=kurtosis['smoothing_sequences']
                                           [k_parms['smooth']],
                energy_window=values['nrg_win'],
                response_file=params['responses'][values['resp']],
                use_polarity=values['polar'],
                n_follow=values['n_follow'])

        gw = params['global_window']
        SNR = params['SNR']
        assoc = params['association']
        val = cls(
            gw_frequency_band=gw['frequency_band'],
            gw_sliding_length=gw['sliding_length'],
            gw_offsets=gw['offsets'],
            SNR_signal_window=SNR['window_lengths']['signal'],
            SNR_noise_window=SNR['window_lengths']['noise'],
            SNR_min_threshold_crossings=SNR['min_threshold_crossings'],
            SNR_quality_thresholds=SNR['pick_quality_thresholds'],
            assoc_cluster_windows=assoc['cluster_windows'],
            station_parameters=sp,
            channel_mapping_rules=ChannelMappingRules())
        if 'dip_rect_thresholds' in params:
            val.dip_rect_thresholds=params['dip_rect_thresholds']
        val.assoc_min_picks = float(assoc.get('min_picks', val.assoc_min_picks))
        val.assoc_max_std = float(assoc.get('max_std', val.assoc_max_std))
        val.assoc_max_std_PtoS = float(assoc.get('max_std_PtoS',
                                       val.assoc_max_std_PtoS))
        if 'responsefiletype' in params:
            val.responsefiletype = params['responsefiletype']
        if 'threshold_parameter' in SNR:
            val.SNR_threshold_parameter = float(SNR['threshold_parameter'])
        if 'distri_secs' in gw:
            val.gw_distri_secs = float(gw['distri_secs'])
        if 'end_cutoff' in gw:
            val.gw_end_cutoff = float(gw['end_cutoff'])
        if 'n_extrema' in gw:
            val.gw_n_extrema = int(gw['n_extrema'])
        if 'extrema_samples' in gw:
            val.gw_extrema_samples = float(gw['extrema_samples'])
        return val


if __name__ == '__main__':
    pass
