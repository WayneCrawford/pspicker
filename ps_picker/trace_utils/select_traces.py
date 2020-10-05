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
    :returns: channel_mappings: sorted dict of ChannelMap with key=station
    """
    channel_map = dict()
    stations = list(set([t.stats.station for t in stream]))
    for station in sorted(stations):
        channel_map[station] = ChannelMap(
            Z=findChannel(stream, station, map_rules.compZ, map_rules.band_order),
            N=findChannel(stream, station, map_rules.compN, map_rules.band_order),
            E=findChannel(stream, station, map_rules.compE, map_rules.band_order),
            H=findChannel(stream, station, map_rules.compH, map_rules.band_order),
            P_write_cmp=map_rules.P_write_cmp,
            S_write_cmp=map_rules.S_write_cmp,
            P_write_phase=map_rules.P_write_phase,
            S_write_phase=map_rules.S_write_phase)

    if verbose:
        print('Channel Map:')
        print(f'{"Station":8s}|'
              + channel_map[stations[0]].__str__(format='table_header'))
        for station in sorted(stations):
            print(f'{station:8s}|'
                  + channel_map[station].__str__(format='table_row'))
    return channel_map


def findChannel(stream, station, cmp_chars, band_order):
    """
    Return the seed ID matching the given station and channel regex

    :param station: station name
    :param cmp_chars: possible characters for component code
    :band order: string listing band codes in order of preference (used only
        if there is more than one trace with valid component code)
    :returns: None if none found, error if more than one found
    """
    id_match = [tr.get_id() for tr in stream.select(station=station)
                if re.search('[' + cmp_chars + ']', tr.stats.channel[-1]) is not None]
    if len(id_match) > 1:
        for band_code in band_order:
            if band_code in [id.split('.')[-1][0] for id in id_match]:
                id_match = [id for id in id_match
                            if id.split['.'][0] == band_code]
    if len(id_match) == 0:
        return None
    assert len(id_match) == 1,\
        'more than one match found for station = {}, cmp = {}'.format(
            station, cmp_chars)
    return id_match[0]


if __name__ == '__main__':
    pass
