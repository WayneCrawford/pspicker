# Generated with SMOP  0.41-beta
# from smop.libsmop import *
import warnings

from obspy.core import Stream, Trace


class PickerLoopParameters():
    """
    """
    def __init__(self,
                 station=None,
                 station_params=None,
                 channel_map=None,
                 datP=None, datS=None, dat_noH=None):
        """
        Parameters associated with work on a single station for one event

        Normally none of the params is provided on initialisation
        :param station: station name
        :param station_params: StationParameters object
        :param channel_map: ChannelMapping object
        """
        self.station = station
        self.station_params = station_params
        self.channel_map = channel_map

    def add_component_traces(self, stream):
        """
        Adds datP, datS, dat_noH, dmin, dmax and t_begin attributes
        """
        self.datP = Stream()
        self.datS = Stream()
        self.dat_noH = Stream()
        self.dmin = 1.e99
        self.dmax = -1.e99
        for comp in self.station_params.P_comp:
            # self.datP += stream.select(id=self.channel_map.__dict__[comp])
            st = stream.select(station=self.station, component=comp)
            if not len(st) == 1:
                txt = f'station={self.station}, comp={comp}'
                if len(st) == 0:
                    warnings.warn(f'No traces found for {txt}')
                    continue
                else:
                    warnings.warn('{:d)} traces found for {}, using first'.
                                  format(len(st), txt))
            self.datP += st[0]
            if self.datP[-1].data.max() > self.dmax:
                self.dmax = self.datP[-1].data.max()
            if self.datP[-1].data.min() < self.dmin:
                self.dmin = self.datP[-1].data.min()
        for comp in self.station_params.S_comp:
            st = stream.select(station=self.station, component=comp)
            if not len(st) == 1:
                txt = f'station={self.station}, comp={comp}'
                if len(st) == 0:
                    warnings.warn(f'No traces found for {txt}')
                    continue
                else:
                    warnings.warn('{:d)} traces found for {}, using first'.
                                  format(len(st), txt))
            self.datS += st[0]
            # self.datP += stream.select(id=self.channel_map.__dict__[comp])
        for tr in stream.select(station=self.station):
            if not tr.stats.channel[-1] == 'H':
                self.dat_noH += tr
        self.t_begin = self.datP[0].stats.starttime

    def index_to_time(self, index, trace=None):
        """
        Returns corresponding to given index

        index: the index
        trace: the trace from which to get the start time and sampling_rate.
            If None, use self.datP[0]
        """
        if trace is None:
            trace = self.datP[0]
        assert isinstance(trace, Trace), f'trace is {type(trace)}, not Trace'
        sr = trace.stats.sampling_rate
        t_begin = trace.stats.starttime
        return t_begin + index/sr


if __name__ == '__main__':
    pass
