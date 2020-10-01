import numpy as np

from obspy.signal.util import smooth as obspy_smooth


def snr_function(M, window_before, window_after):
    """
    Return estimated signal to noise ratio (dB) at time t

    :param M: energy (obspy Trace object)
    :param window_before: length of "noise" window (seconds before a
        given time)
    :param window_after: length of "signal" window (seconds after a
        given time)
    :returns: 20*log10(M(t...t+window_after)/M(t-window_before,t))
    """
    sr = M.stats.sampling_rate
    samps_beforewind = round(sr * window_before)
    samps_afterwind = round(sr * window_after)
    w_before = M.copy()
    w_after = M.copy()
    w_before.data = obspy_smooth(M.data, samps_beforewind)
    w_after.data = obspy_smooth(M.data, samps_afterwind)
    # Shift times correctly from half the size of both windows
    w_before.stats.starttime += window_before/2
    w_after.stats.starttime -= window_after/2
    global_start = w_before.stats.starttime
    global_end = w_after.stats.endtime
    w_before.trim(starttime=global_start, endtime=global_end)
    w_after.trim(starttime=global_start, endtime=global_end)
    snr = w_after.copy()
    snr.data /= w_before.data
    snr.data = 20 * np.log10(snr.data)
    # B_before = concat([[Aa(arange(wa,end()),arange())],[zeros(wa - 1,n)]])
    # B_after = concat([[Aa(arange(wa,end()),arange())],[zeros(wa - 1,n)]])
    # Bb=copy(Ab)
    # snr=abs(Ba / Bb)
    # snr=dot(20,log10(snr))
    # snr[arange(1,wb),arange()]=0
    return snr


if __name__ == '__main__':
    pass
