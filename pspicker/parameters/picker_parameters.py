# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# readmain.m
import yaml

from .channel_mapping_rules import ChannelMappingRules
from .station_parameters import StationParameters
from .SNR_parameters import SNRParameters
from .polarity_parameters import PolarityParameters
from .associator_parameters import AssociatorParameters
from .global_window_parameters import GlobalWindowParameters
from ..timer import Timer


class PickerParameters():
    """
    Picker Parameters
    """
    def __init__(self, global_window, SNR, polarity={}, channel_parameters={},
                 association={}, response_file_type='', station_parameters={},
                 stations={}):
        """
        Initialize Picker Parameters

        :param gw: GlobalWindowParameters dictionary source
        :param SNR: SNRParameters dictionary source
        :param polarity: PolarityParameters dict source
        :param assoc: AssociatorParameters dict source
        :param station_parameters: dict with key=station type,
                                   items=StationParameters object
        :param stations: dict with key=station name,
                                   items={station type, resp_file_name}
        :param channel_parameters: ChannelMappingRules dictionary source
        :param response_file_type: Format of the response file(s).  ['GSE',
            'JSON_PZ', 'SACPZ', 'STATIONXML' or
            empty.  If empty use Baillardesque PoleZero format.
        """
        self.gw = GlobalWindowParameters(**global_window)
        self.SNR = SNRParameters(**SNR)
        self.polarity = PolarityParameters(**polarity)
        self.channel_mapping_rules = ChannelMappingRules(**channel_parameters)
        self.assoc = AssociatorParameters(**association)
        self.response_file_type = response_file_type
        
        self.station_parameters = {}
        for station, values in stations.items():
                temp = station_parameters[values['parameters']]
                temp['resp_file'] = values['resp_file']
                self.station_parameters[station] = StationParameters.from_dict(temp)

    @property
    def stations(self):
        return [s for s in self.station_parameters.keys()]

    @property
    def n_stations(self):
        return len(self.station_parameters)

    def __str__(self):
        str = "PickerParameters:\n"
        str += f"    gw = {self.gw}\n"
        str += f"    SNR = {self.SNR}\n"
        str += f"    polarity = {self.polarity}\n"
        str += f"    assoc = {self.assoc}\n"
        str += f"    stations = {','.join(self.stations)}\n"
        str += f"    response_file_type = '{self.response_file_type}'\n"
        str += f"    channel_mapping_rules = {self.channel_mapping_rules}\n"
        return str

    @classmethod
    def from_yaml_file(cls, filename):
        """
        Read in parameters from a YAML file
        """
        with Timer(text="Read picker parameter file: {:0.4f}s"):
            with open(filename, 'r') as fic:
                params = yaml.load(fic)

            return cls(**params)


if __name__ == '__main__':
    pass

#     def __init__(self,
#                  global_window,
#                  SNR,
#                  assoc,
#                  station_parameters,
#                  channel_mapping_rules,
#                  polarity,
#                  response_file_type=''):
#         """
#         Initialize Picker Parameters
# 
#         :param gw: GlobalWindowParameters object (control inital selection
#             of pick window)
#         :param SNR: SNRParameters object
#         :param polarity: PolarityParameters object
#         :param assoc: AssociatorParameters object
#         :param station_parameters: dict with key=station names,
#             items=StationParameters object
#         :param channel_mapping_rules: ChannelMappingRules object
#         :param response_file_type: Format of the response file(s).  'GSE' or
#             empty.  If empty use Baillardesque PoleZero format.
#         """
#         self.gw = gw
#         self.SNR = SNR
#         self.polarity = polarity
#         self.assoc = assoc
#         self.station_parameters = station_parameters
#         self.channel_mapping_rules = channel_mapping_rules
#         self.response_file_type = response_file_type
#         print(self)

#     @classmethod
#     def from_yaml_file(cls, filename):
#         """
#         Read in parameters from a YAML file
#         """
#         with Timer(text="Read picker parameter file: {:0.4f}s"):
#             with open(filename, 'r') as fic:
#                 params = yaml.load(fic)
# 
#             return cls(**params)
#             # Fill in station parameters
#             # with Timer(text="    Station Parameters: {:0.4f}s"):
#             sp = dict()
#             for station, values in params['stations'].items():
#                 temp = params['station_parameters'][values['parameters']]
#                 temp['resp_file'] = values['resp_file']
#                 sp[station] = StationParameters.from_dict(temp)
#             cmr = ChannelMappingRules.from_dict(params.get(
#                 'channel_parameters', {}))
#             # with Timer(text="    Global Window Parameters: {:0.4f}s"):
#             gw = GlobalWindowParameters.from_dict(params['global_window'])
#             # with Timer(text="    SNR Parameters: {:0.4f}s"):
#             SNR = SNRParameters.from_dict(params['SNR'])
#             # with Timer(text="    Polarity Parameters: {:0.4f}s"):
#             polarity = PolarityParameters.from_dict(params.get('polarity', {}))
#             # with Timer(text="    Associator Parameters: {:0.4f}s"):
#             assoc = AssociatorParameters.from_dict(params['association'])
#             # with Timer(text="    Overall Parameters: {:0.4f}s"):
#             val = cls(
#                 gw=gw,
#                 SNR=SNR,
#                 assoc=assoc,
#                 polarity=polarity,
#                 station_parameters=sp,
#                 channel_mapping_rules=cmr)
#             if 'response_file_type' in params:
#                 val.response_file_type = params['response_file_type']
#         return val

