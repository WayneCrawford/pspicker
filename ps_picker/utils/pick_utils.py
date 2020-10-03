"""
Routines involving lists of obspy Picks
"""


def picks_ps_times(picks):
    """
    Return ps_times and p_time for events having P and S picks

    :param picks: list of obspy Pick objects
    :returns: delays, P_times
    :rtype: list, list of UTCDateTime
    """
    matches = picks_matched_stations(picks)
    if len(matches) == 0:
        return None, None
    return ([(x['pickS'].time - x['pickP'].time) for x in matches],
            [x['pickP'].time for x in matches])


def picks_matched_stations(picks):
    """
    Return list of stations with P AND S picks

    :param picks: list of obspy Pick objects
    :returns: list of {'station', 'pickP', 'pickS'}
    """
    picksP = [p for p in picks if p.phase_hint[0] == 'P']
    picksS = [p for p in picks if p.phase_hint[0] == 'S']
    matches = []
    for pickP in picksP:
        for pickS in picksS:
            if pickP.station == pickS.station:
                matches.append({'station': pickP.station,
                                'pickP': pickP, 'pickS': pickS})
                break
    return matches

# def _old_match_pick_stations(picksP, picksS):
#     """
#     Return list of stations with P AND S picks, and the indices of the picks
#
#     :param picksP: list of P picks
#     :param picksS: list of S picks
#     :returns: list of {'station', 'iP', 'iS'}
#     """
#     matches = []
#     iP = 0
#     for pickP in picksP:
#         iS = 0
#         for pickS in picksS:
#             if pickP.station == pickS.station:
#                 matches.append({'station': pickP.station,
#                                 'iP': iP, 'iS': iS})
#                 break
#             iS += 1
#         iP += 1
#     return matches
