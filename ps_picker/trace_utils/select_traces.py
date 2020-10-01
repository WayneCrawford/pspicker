# from readMSEEDTraces.m  ?
import re

from ..utils import ChannelMap


def select_traces(stream, map_rules, verbose=1):
    """
    select traces and other stuff

    based on old readMSEEDTraces()

    :param stream: obspy stream
    :param map_rules: ChannelMappingRules object
    :param verbose: talk
    :returns: channel_mappings: dict with key=station, subkeys=Z,E,N,H,
        P_write, S_write
    """
    channel_map = dict()
    stations = list(set([t.stats.station for t in stream]))
    for station in stations:
        channel_map[station] = ChannelMap(
            Z=findChannel(stream, station, map_rules.compZ),
            N=findChannel(stream, station, map_rules.compN),
            E=findChannel(stream, station, map_rules.compE),
            H=findChannel(stream, station, map_rules.compH),
            P_write_to=map_rules.putPick_P_Comp,
            S_write_to=map_rules.putPick_S_Comp)

    if verbose:
        print('Channel Map:')
        fmt = '{:8s}|{:^16s}|{:^16s}|{:^16s}|{:^16s}|'
        print(fmt.format('Station', 'Z', 'N', 'E', 'H'))
        print('{:-<8s}+{:-<16s}+{:-<16s}+{:-<16s}+{:-<16s}|'.format(
            '', '', '', '', ''))
        for station in sorted(stations):
            H = channel_map[station].H
            if H is None:
                H = 'None'
            print(fmt.format(station,
                             channel_map[station].Z,
                             channel_map[station].N,
                             channel_map[station].E,
                             H))
    return channel_map


def findChannel(stream, station, compRegEx):
    """
    Return the seed ID matching the given station and channel regex

    :returns: None if none found, error if more than one found
    """
    id_match = [tr.get_id() for tr in stream.select(station=station)
                if re.search(compRegEx, tr.stats.channel) is not None]
    if len(id_match) == 0:
        return None
    assert len(id_match) == 1,\
        'more than one match found for station = {}, channel = {}'.format(
            station, compRegEx)
    return id_match[0]


if __name__ == '__main__':
    pass
