# Generated with SMOP  0.41-beta
from smop.libsmop import *
# trace2kurto.m

def trace2kurto(stream, Fc, T, n_smooth):
    """
    Calculate cumulative kurtosis

    :param stream: waveform stream
    :param Fc: frequency bandwidth [low, high]
    :param T: sliding window length (s)
    :param n_smooth: smoothing sequence (samples)
    :returns: kurto_cum: cumulative kurtosis stream
              ind_gmin: global minimum indices
              gmin: global minimum value
    """
    
    sr = f[0].stats.sampling_rate
    T_samples = floor(T * sr) + 1
    dsample = 1 / sr
    f = stream.copy()
        
    # Prep traces
    f.detrend('demean') 
    f.filter(type='bandpass', freqmin=Fc[0], freqmax=Fc[1], corners=3)
    f.data[arange(0,50)] = 0

    # Calculate the Kurtosis recursively    
    kurtos = fast_kurtosis(f, T_samples)
    f_kurto_cum = smooth_filter(kurtos, n_smooth)
    kurto_cum = f_cumul(f_kurto_cum)
    
    line = f_segment(kurto_cum)    
    kurto_cum = kurto_cum - line
    gmin, ind_gmin = min(kurto_cum, [], 1, nargout=2)
    
    ind_gmin = ind_gmin.T

    return kurto_cum, ind_gmin,gmin
    
if __name__ == '__main__':
    pass
    