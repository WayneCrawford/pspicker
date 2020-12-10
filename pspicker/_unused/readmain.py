# Generated with SMOP  0.41-beta
from smop.libsmop import *
import yaml
import os.path
# readmain.m

    # Function that outputs reapaths and wavpaths between two dates 
# and gives the directory where the new reafiles should be written
# NB: Do an AUTOREG (Seisan) before launching picker to put data in the
# database (REA+WAV)
    
class StationParameters():
    """
    Station Parameters
    """
    def __init__(self, P_comp, S_comp, f_energy, frequencies, window_lengths,
                 smoothing_sequence, energy_window, response, polarity=False,
                 lag=999, n_follow=2)
        """
        P_Comp: String of components used to search for P arrival & amplitude
        S_comp: String of components used to search for S arrival & amplitude
        f_energy:  frequency range in which to search for SNR [lo, high]
        :param frequencies: kurtosis frequency bands (list of [lo, high]s)
        :param window_lengths: kutosis window lengths (list of seconds)
        :param smoothing_sequence: kurtosis smooting sequence (increasing
            list of samples(?))
        :param energy_window: only look at data from t-winlen to t, where t
            is the time of the peak waveform energy.  If == 0, don't use energy
            criteria.
        :param response: response file name
        :param polarity: Use polarity analysis to identify phases?
        :param lag: ??? (unused?)
        :param n_follow: number of extrema to follow (1 or 2).  Generally use
            2 (S and P) unless data are problematic (some OBS data, for example)
        """
        self.P_comp = P_comp
        self.S_comp = S_comp
        self.f_energy = f_energy
        self.frequencies = frequencies
        self.window_lengths = window_lengths
        self.smoothing_sequence = smoothing_sequence
        self.energy_window = energy_window
        self.response = response
        self.polarity = polarity
        self.lag = lag
        self.n_follow = n_follow
    

class PickerParameters():
    """
    Picker Parameters
    """
    def __init__(self, rewindow_frequency, rewindow_sliding_length,
                 rewindow_offsets, SNR_window_offsets, SNR_quality_threshold,
                 dip_rect_thresholds, cluster_windows, SNR_quality_thresholds,
                 amplitude_remove, kurtosis_frequency_bands, kurtosis_window_lengths,
                 kurtosis_smoothing_sequences, responses, station_parameters,
                 base_name, start_year, start_month, end_year, end_month,
                 input_path, output_path='./Sfile_directory'):
        """
        Initialize Picker Parameters
        
        :param rewindow_frequency: The frequency range for rewindowing [low, high]
        :param rewindow_sliding_length: Length in seconds of sliding window for
            Kurtosis calculation
        :param rewindow_offsets: Offsets in seconds for global (all station based)
            rewindowing [left, right]
        :param SNR_window_offsets: Offsets in seconds for SNR calculation
            [left, right]
        :param SNR_quality_threshold: SNR trace quality threshold
        :param dip_rect_thresholds: DipRect thresholds {'P': value, 'S': value}
        :param cluster_windows: Windows in seconds for clustering rejection
            {'P': value, 'S':value}
        :param pick_quality_thresholds: SNR threshold for pick quality assignment
            quality = {'3':value, '2':value, '1':value, '0':value]
            name changed from SNR quality thresholds
        :param amplitude_remove: Remove old Amplitudes in Amplitude picking mode?
            bool, ignored?
        :param station_parameters: dict with key=station names,
            items=StationParameters object

        Note: The base_name, input_path and start/end dates are intended to
              allow processing of a big chunk of data, automatically determining
              the database and waveform directories, but PSPicker currently just
              works on one event, which is provided in its argument lis.
        :param base_name: SEISAN database name
        :param input_path: Directory in which the SEISAN /WAV directory is found
        :param output_path: Output directory
        :param start_year: start year of data to process
        :param end_year: end year of data to process
        :param start_month: start month of data to process
        :param end_month: end month of data to process
        """
        self.output_path = output_path
        self.rewindow_frequency = rewindow_frequency
        self.rewindow_sliding_length = rewindow_sliding_length
        self.rewindow_offsets = rewindow_offsets
        self.SNR_window_offsets = SNR_window_offsets
        self.SNR_quality_threshold = SNR_quality_threshold
        self.dip_rect_thresholds = dip_rect_thresholds
        self.cluster_windows = cluster_windows
        self.pick_quality_thresholds = pick_quality_thresholds
        self.amplitude_remove = amplitude_remove
        self.station_parameters = station_parameters
        # Unused?
        self.base_name = base_name
        self.input_path = input_path
        self.start_year = start_year
        self.start_month = start_month
        self.end_year = end_year
        self.end_month = end_month
        self.rea_paths = _make_seisan_paths('REA')
        self.wav_paths = _make_seisan_paths('WAV')
        self.n_months = end_month - start_month + 1
        if self.end_year != self.start_year:
        self.n_months += 12*(self.end_year-self.start_year)

        
    def _make_seisan_paths(self,name)
        """
        make seisan paths covering start_date to end_date
        :param name: name of the subdirectory ('REA' or 'WAV')
        """
        paths = []
        year = self.start_year
        month = self.start_month
        while year <= self.end_year:
            if year == self.end_year and month > self.end_month:
                break
            if month > 12:
                year += 1
                month = 1
                continue
            paths.extend(os.path.join(self.input_path, 'REA', self.base_name,
                                      f'{year:04d}', f'{month:02d}', ''))
            month += 1
        return paths
    
    @classdef
    from_yaml_file(cls, filename):
    """
    Read in parameters from a YAML file
    """
    with open(filename, 'r') as fic:
        params = yaml.load(fic)

    # Fill in station parameters
    sp = dict()
    for station, values in sp.items():
        k_parms = values['k_parms']
        sp[station] = StationParameters(
            P_comp=values['P_comp'],
            S_comp=values['S_comp'],
            f_energy=values['f_nrg'],
            frequencies=params['kurtosis']['frequency_bands'][k_parms['freqs']]
            window_lengths=params['kurtosis']['window_lengths'][k_parms['wind']]
            smoothing_sequence=params['kurtosis']['smoothing_sequences']
                                     [k_parms['smooth']]
            energy_window=values['nrg_win']
            response=params['responses'][values['resp']]
            polarity=values['polar']
            lag=values['lag']
            n_follow=values['n_follow'])

    start_year = int(params['dates'][0][:4])
    start_month = int(params['dates'][0][4:])
    end_year = int(params['dates'][1][:4])
    end_month = int(params['dates'][1][4:])
    
    val = cls(rewindow_frequency=params['rewindowing']['frequency'],
              rewindow_sliding_length=params['rewindowing']['sliding_length'],
              rewindow_offsets=params['rewindowing']['offsets'],
              SNR_window_offsets=params['SNR']['window_offsets'],
              SNR_quality_threshold=params['SNR']['quality_threshold'],
              pick_quality_thresholds=params['SNR']['pick_quality_thresholds'],
              dip_rect_thresholds=params['dip_rect_thresholds'],
              cluster_windows=params['cluster_windows'],
              amplitude_remove=params['amplitude_remove'],
              station_parameters=sp,
              base_name=params['base_name'],
              start_year=start_year,
              start_month=start_month,
              end_year=end_year,
              end_month=end_month,
              input_path=params['input_path'],
              output_path=params['output_path'])
    return val

    # Edit all SEISAN REA and WAV paths were data are

@function
def readmain(filename=None,*args,**kwargs):
    # varargin = readmain.varargin
    # nargin = readmain.nargin

    # READMAIN reads PSPicker parameter file
    
    # Usage: PickerParam=readmain(filename)
    
    # The parameter file is divided into sections started by a comment
    #   line (first character '#') and ending with an empty line
    
    # The sections are (order is important!):
    #   1:  The SEISAN database name
    #   2:  The directory in which the SEISAN /WAV directory is found
    #   3:  The start and endtime (YYYYMM) (returns appropriate REA and WAV
    #        directories, but PSPicker currently just processes one file,
    #        whose name and path are specified in its argument list).
    #   4:  The output directory (ignored, PSPicker uses './Sfile_directory')
    #   5:  The frequency range for rewindowing [low high]
    #   6:  Length of sliding window for Kurtosis calculation (seconds)
    #   7:  Left and Right offsets for global (all station based) rewindowing (sec)
    #   8:  Left and Right window for SNR calculation
    #   9:  SNR trace quality threshold
    #   10: DipRect P and S thresholds
    #   11: P and S windows for clustering rejection (sec)
    #   12: Pick SNR threshold for pick quality assignment (3 2 1 0)
    #   13: Remove old Amplitudes in Amplitude picking mode? (0 or 1, ignored?)
    #   14: Frequency Bands (Hz) headed by a character code
    #   15: Window Lengths (s) headed by a character code
    #   16: Smoothings (samples) headed by a character code
    #   17: Response Files' names, headed by a character code
    #   18: Station-dependent parameters
    #       The station dependent parameter liness consist of the following columns
    #       1:  Station name
    #       2:  Components used to search for P arrival & amplitude
    #       3:  Components used to search for S arrival & amplitude
    #       4:  f_energy
    #       5:  frequencies (a code from the "Frequency Bands" section)
    #       6:  window      (code from the "Window Lengths" section)
    #       7:  smoothing   (code from the "Smoothings" section)
    #       8:  polarity (1 or 0: use polarity analysis to identify phases)
    #       9:  lag (???, unused)
    #       10: energy_winlen: only look at data from `t-winlen` to `t`,
    #           where `winlen1 is in seconds and `t` is the time of the 
    #           peak waveform energy (`winlen`= 0 means don't use energy criteria)
    #       11: n_follow: number of extrema to follow (1 or 2): generally
    #                   use 2 (S&P) unless data are problematic (some OBS
    #                   data, for example).
    #       12: response    (code from the "Response Files" section)
    #
    #   Note: The first 3 sections (SEISAN database name, SEISAN directory
    #       and start and end dates) are intended to allow processing of
    #       a big chunk of data, automatically determining the database
    #       and waveform directories, but PSPicker currently just works
    #       on one event, which is provided in its argument lis.

    with open(filename, 'r') as fic:
        params = yaml.load(fic)

    key_list = ['base', 'input_path', 'dates', 'output_path',
                'rewindow_frequency', 'rewindow_sliding_length',
                'rewindow_offsets', 'SNR_window_offsets',
                'SNR_quality_threshold', 'dip_rect_thresholds',
                'cluster_windows', 'SNR_quality_thresholds',
                'amplitude_remove', 'kurtosis_frequency_bands',
                'kurtosis_window_lengths', 'kurtosis_smoothing_sequences',
                'responses', 'station_parameters']
    for key in key_list:
        assert key in params, f'no "{key}" key in "{filename}"'

    # Fill in station parameters
    sp = params['station_parameters']
    for station in sp:
        k_parms = sp[station]['k_parms']
        sp[station]['k_freqs'] = params['kurtosis_frequency_bands'][k_parms['freqs']]
        sp[station]['k_wind'] = params['kurtosis_window_lengths'][k_parms['wind']]
        sp[station]['k_smooth'] = params['kurtosis_smoothing_sequences'][k_parms['smooth']]
        sp[station]['resp'] = params['responses'][sp[station]['resp']]

    # Get total number of months to be processed
    start_year = int(params['dates'][0][:4])
    start_month = int(params['dates'][0][4:])
    end_year = int(params['dates'][1][:4])
    end_month = int(params['dates'][1][4:])

    n_months = end_month - start_month + 1
    if end_year != start_year:
        n_months += 12*(end_year-start_year)

    # Edit all SEISAN REA and WAV paths were data are

    params['rea_paths'] = []
    params['wav_paths'] = []
    year = start_year
    month = start_month
    while year <= end_year:
        if year == end_year and month > end_month:
            break
        if month > 12:
            year += 1
            month = 1
            continue
        params['rea_paths'].extend(os.path.join(params['input_path'], 'REA',
                                                params['base'], f'{year:04d}',
                                                f'{month:02d}', ''))
        params['wav_paths'].extend(os.path.join(params['input_path'], 'WAV',
                                                params['base'], f'{year:04d}',
                                                f'{month:02d}', ''))
        month += 1

    return params


if __name__ == '__main__':
    pass


# The old reading and assigning section
#     Base=C[1][1]
#     Input_path=C[2][1]
#     Dates=textscan(C[3][1],'%4.0f%2.0f %4.0f%2.0f')
#     Start_year=Dates[1]
#     Start_month=Dates[2]
#     End_year=Dates[3]
#     End_month=Dates[4]
#     Output_path=C[4][1]
#     R_freq=str2num(C[5][1])
#     R_slid=str2num(C[6][1])
#     R_wind=textscan(C[7][1],'%f %f')
#     R_wind=concat([R_wind[arange()]])
#     SNR_wind=textscan(C[8][1],'%f %f')
#     SNR_wind=concat([SNR_wind[arange()]])
#     SNR_quality=str2num(C[9][1])
#     DipRect=textscan(C[10][1],'%f %f')
#     DipRect=concat([DipRect[arange()]])
#     PS_clust=textscan(C[11][1],'%f %f')
#     PS_clust=concat([PS_clust[arange()]])
#     SNR_thres=str2num(C[12][1])
#     Amplitude_remove=str2num(C[13][1])
#     Frequencies=C[14]
#     Windows=C[15]
#     Smoothings=C[16]
#     Amplitude=C[17]
#     Stat_param=C[18]
#     # Assign station dependent parameters to variables

#     for i in arange(1,length(Stat_param)).reshape(-1):
#         clear('scratch')
#         line=Stat_param[i]
# # readmain.m:121
#         left_br=regexpi(line,'[')
# # readmain.m:122
#         right_br=regexpi(line,']')
# # readmain.m:123
#         F_energy=str2num(line(arange(left_br(1),right_br(1))))
# # readmain.m:124
#         line[arange(left_br(1),right_br(1))]=[]
# # readmain.m:125
#         scratch=textscan(line,'%s %s %s %s %s %s %f %f %f %f %s')
# # readmain.m:126
#         Station=scratch[1]
# # readmain.m:127
#         P_comp=scratch[2]
# # readmain.m:128
#         S_comp=scratch[3]
# # readmain.m:129
#         Kurto_F=scratch[4]
# # readmain.m:130
#         Kurto_W=scratch[5]
# # readmain.m:131
#         Kurto_S=scratch[6]
# # readmain.m:132
#         Polarity_flag=scratch[7]
# # readmain.m:133
#         Lag=scratch[8]
# # readmain.m:134
#         Energy_recond=scratch[9]
# # readmain.m:135
#         N_follow=scratch[10]
# # readmain.m:136
#         Ampli=scratch[11]
# # readmain.m:137
#         ind_F=strcmpi(Kurto_F,cellfun(lambda x=None: x(1),Frequencies,'un',0))
# # readmain.m:139
#         Fref=str2num(Frequencies[ind_F](arange(2,end())))
# # readmain.m:140
#         ind_W=strcmpi(Kurto_W,cellfun(lambda x=None: x(1),Windows,'un',0))
# # readmain.m:141
#         Winf=str2num(Windows[ind_W](arange(2,end())))
# # readmain.m:142
#         ind_S=strcmpi(Kurto_S,cellfun(lambda x=None: x(1),Smoothings,'un',0))
# # readmain.m:143
#         Smoof=str2num(Smoothings[ind_S](arange(2,end())))
# # readmain.m:144
#         ind_A=strcmpi(Ampli,cellfun(lambda x=None: x(1),Amplitude,'un',0))
# # readmain.m:145
#         Ampf=strtrim(Amplitude[ind_A](arange(2,end())))
# # readmain.m:146
#         Station_param[i]=struct('Station',Station,'P_comp',P_comp,'S_comp',S_comp,'F_energy',F_energy,'Kurto_F',Fref,'Kurto_W',Winf,'Kurto_S',Smoof,'Polarity_flag',Polarity_flag,'Lag',Lag,'Energy_recond',Energy_recond,'N_follow',N_follow,'Amplitude',Ampf)
# # readmain.m:148
