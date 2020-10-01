from scipy.signal import find_peaks


def pk2pk(stream, pick_time, before_pick, after_pick):
    """
    Returns the time, period and amplitude of the maximum peak-peak

    :param stream: waveform data
    :param pick_time: time of reference (usually a pick)
    :param before_pick = seconds before pick for left end of window
    :param after_pick = seconds after pick for right end of window
    :returns: {'period': value, 'amplitude': value, 'time': UTCDateTime]
    """
    Amp = {'time': None, 'amplitude': 0, 'period': None}
    start_time = pick_time - before_pick
    end_time = pick_time + after_pick
    for tr in stream:
        window = tr.copy().trim(start_time, end_time).data
        sr = tr.stats.sampling_rate
        i_mins = find_peaks(-1 * window.data)
        i_maxs = find_peaks(window.data)

        # Make sure there is an i_max after the last i_min
        if i_mins[-1] == len(i_mins)-1:
            i_mins.pop(-1)
        if i_maxs[-1] <= i_mins[-1]:
            i_maxs.extend(len(i_mins)-1)

        for imin in i_mins:
            imax = i_maxs(i_maxs > imin)[0]
            amplitude = window.data[imax] - window.data[imin]
            if amplitude > Amp['amplitude']:
                Amp = {'time': start_time + (imin / sr),
                       'amplitude': amplitude,
                       'period': (imax - imin) / sr}
    return Amp


if __name__ == '__main__':
    pass
