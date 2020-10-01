import numpy as np


class PickerRunParameters():
    """
    Parameters associated with the run of one event
    """
    def __init__(self,
                 rea_name,
                 wavefile=None,
                 stream=None,
                 channel_maps=None,
                 overall_distri=np.array(0),
                 global_first_time=None,
                 global_last_time=None,
                 plotter=None):
        """
        :param rea_name: database file to read
        :param wavefile: name of the file containing the waveforms
        :param stream: all of the data waveforms
        :param channel_maps: a dict of ChannelMapping objects with station
            names as key
        :param overall_distri: numpy array containing offset of each extrema
        :param t_begin: time of reference for overall_distri
        :param global_first_time: never look before this time
        :param global_last_time: never look after this time
        :param plotter: PSPickerPlotter object for making plots
        """
        self.rea_name = rea_name
        self.channel_maps = channel_maps
        self.wavefile = wavefile
        self.stream = stream
        self.overall_distri = overall_distri
        self.global_first_time = global_first_time
        self.global_last_time = global_last_time
        self.plotter = plotter

    @property
    def stations(self):
        return [s for s in self.channel_maps.keys()]

    @property
    def n_stations(self):
        return len(self.stations)


if __name__ == '__main__':
    pass
