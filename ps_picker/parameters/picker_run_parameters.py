class PickerRunParameters():
    """
    Parameters associated with the run of one event
    """
    def __init__(self,
                 database_filename,
                 wavefile=None,
                 stream=None,
                 channel_maps=None,
                 first_time=None,
                 last_time=None,
                 t_begin=None):
        """
        :param database_filename: database file to read
        :param wavefile: name of the file containing the waveforms
        :param stream: all of the data waveforms
        :param channel_maps: a dict of ChannelMapping objects with station
            names as key
        :param t_begin: time of reference for overall_distri
        :param first_time: never look before this time
        :param last_time: never look after this time
        """
        self.database_filename = database_filename
        self.channel_maps = channel_maps
        self.wavefile = wavefile
        self.stream = stream
        self.first_time = first_time
        self.last_time = last_time
        self.t_begin = t_begin

    @property
    def stations(self):
        return [s for s in self.channel_maps.keys()]

    @property
    def n_stations(self):
        return len(self.stations)


if __name__ == '__main__':
    pass
