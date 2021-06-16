# from readMSEEDTraces.m  ?
from ..channel_map import ChannelMap
from ..logger import log


def select_traces(stream, map_rules, debug=False):
    """
    Map traces and phases

    selects traces corresponding to "Z", "N", "E" and "H", names the
    components to write picks to and defines the phase names
    :param stream: obspy stream
    :param map_rules: ChannelMappingRules object
    :param verbose: talk
    :returns: channel_mappings: sorted dict of ChannelMap with key=station
    """
    channel_map = dict()
    stations = list(set([t.stats.station for t in stream]))
    mr = map_rules
    for station in sorted(stations):
        Z = mr.seed_id('Z',stream, station)
        if Z is None:
            sta_st = stream.select(station=station)
            log('No Z channel found for station {}: {} not in "{}"'.format(
                station, [t.stats.channel[-1] for t in sta_st], mr.component_orientation_codes.Z),
                'error')
        else:
            channel_map[station] = ChannelMap(
                Z=mr.seed_id('Z',stream, station),
                N=mr.seed_id('N',stream, station),
                E=mr.seed_id('E',stream, station),
                H=mr.seed_id('H',stream, station),
                P_write_cmp=mr.write_components_phases.P[0],
                S_write_cmp=mr.write_components_phases.S[0],
                P_write_phase=mr.write_components_phases.P[1],
                S_write_phase=mr.write_components_phases.S[1])

    if debug:
        print('Channel Map:')
        print(f'{"Station":8s}|'
              + channel_map[stations[0]].__str__(format='table_header'))
        for station in sorted(stations):
            print(f'{station:8s}|'
                  + channel_map[station].__str__(format='table_row'))
    return channel_map


if __name__ == '__main__':
    pass
