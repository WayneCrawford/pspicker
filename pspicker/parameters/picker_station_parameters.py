# Generated with SMOP  0.41-beta
# from smop.libsmop import *
import warnings

from obspy.core import Stream  # , Trace


class PickerStationParameters():
    """
    Parameters associated with work on a single station for one event
    """
    def __init__(self,
                 station=None,
                 station_params=None,
                 channel_map=None,
                 stream=None):
        """
        Initiate loop parameters

        :param station: station name
        :param station_params: StationParameters object
        :param channel_map: ChannelMapping object
        :param stream: obspy Stream containing appropriate traces
        """
        self.station = station
        self.station_params = station_params
        self.channel_map = channel_map
        self.datP = self._get_traces(stream, station_params.picking_components.P)
        self.datS = self._get_traces(stream, station_params.picking_components.S)
        # print(self.datP)
        # print(self.datS)
        self.dat_noH = self._get_noH_traces(stream)
        self.data_limits = self._find_limits()
        self.t_begin = self._find_first_time()

    def _get_traces(self, stream, comp_list):
        """
        Adds datP or datS

        :param stream: Stream of all Traces
        :param comp_list: list of components to include in this stream
        """
        out_stream = Stream()
        # make sure 'Z' is at end of comp_list (for polarity analysis)
        if 'Z' in comp_list and len(comp_list) > 1:
            comp_list = comp_list.replace('Z', '') + 'Z'
        for comp in comp_list:
            st = stream.select(id=getattr(self.channel_map, comp))
            # st = stream.select(station=self.station, component=comp)
            if not len(st) == 1:
                txt = f'station={self.station}, comp={comp}'
                if len(st) == 0:
                    warnings.warn(f'No traces found for {txt}')
                    continue
                else:
                    warnings.warn('{:d)} traces found for {}, using first'.
                                  format(len(st), txt))
            out_stream += st[0]
        return out_stream

    def _get_noH_traces(self, stream):
        """
        Adds dat_noH
        """
        out_stream = Stream()
        for tr in stream.select(station=self.station):
            if not tr.stats.channel[-1] == 'H':
                out_stream += tr
        return out_stream

    def _find_limits(self):
        """
        find limits to self.datP, self.datS and self.dat_noH

        Uses only self.dat_noH because the others should be included
        """
        assert self.dat_noH is not None, 'self.dat_noH is None'
        assert len(self.dat_noH) > 0, 'no traces in self.dat_noH'
        dmin = 1.e99
        dmax = -1.e99
        for tr in self.dat_noH:
            if tr.data.max() > dmax:
                dmax = tr.data.max()
            if tr.data.min() < dmin:
                dmin = tr.data.min()
        return (dmin, dmax)

    def _find_first_time(self):
        """
        Finds first time of all traces

        Uses only self.dat_noH because the others should be included
        """
        assert self.dat_noH is not None, 'self.dat_noH is None'
        assert len(self.dat_noH) > 0, 'no traces in self.dat_noH'

        first_time = self.dat_noH[0].stats.starttime
        for tr in self.dat_noH[1:]:
            if tr.stats.starttime < first_time:
                first_time = tr.stats.starttime
        return first_time

    # def index_to_time(self, index, trace=None):
    #     """
    #     Returns corresponding to given index
    #
    #     index: the index
    #     trace: the trace from which to get the start time and sampling_rate.
    #         If None, use self.datP[0]
    #     """
    #     if trace is None:
    #         trace = self.datP[0]
    #     assert isinstance(trace, Trace), f'trace is {type(trace)}, not Trace'
    #     sr = trace.stats.sampling_rate
    #     t_begin = trace.stats.starttime
    #     return t_begin + index/sr


if __name__ == '__main__':
    pass
