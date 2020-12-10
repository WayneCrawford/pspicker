# Generated with SMOP  0.41-beta
from smop.libsmop import *
# true_amplitude.m

    
@function
def true_amplitude(norm_cons=None,dig_gain=None,cali_gain=None,f_ref=None,poles=None,zes=None,data=None,rsample=None,*args,**kwargs):
    varargin = true_amplitude.varargin
    nargin = true_amplitude.nargin

    # TRUE_AMPLITUDE calculates frequency response function T(f) based on poles
# and zeros
    
    # Inputs:   norm_cons   normalization constant in V/nm
#           dig_gain    digitizer gain in count/V, usually this value is
#                       big due to digitizer senstivity (ex 1e6 count/V)
#           cali_gain   calibration gain in nm/c, this value is related to
#                       a fequence reference value, often 1 Hz.
#           f_ref       frequence reference value
#           poles       expressed in rad/s
#           zes         zeros expressed in rad/s
    
    # For more information see "Of poles and zeros" by Frank Sherbaum
# GSE format are displacement response file, that means we consider the
# output being in counts and the input in displacement (nm)
    
    # The constants(norm_cons and dig_gain) are here to convert adimsensionnal
# FRF to count/nm FRF
# Calibration gain: Suppose the calibration gain is 1/2 nm/c at reference period 1s
# that means for frequence 1 Hz, |T(w_ref)|=2 c/nm=1/gc
    
    # |Tnew(w_ref)|=c*|Told(w_ref)|=1/gc  =>  c=1/[|Told(w_ref)|*gc]
    
    # So the overall FRF is FRF=norm_cons*dig_gain*c*FRF_old
    
    m=size(data,1)
# true_amplitude.m:30
    n=size(data,2)
# true_amplitude.m:31
    NFFT=2 ** nextpow2(m)
# true_amplitude.m:32
    Y=fft(data,NFFT)
# true_amplitude.m:33
    f=dot(rsample / 2,linspace(0,1,NFFT / 2))
# true_amplitude.m:34
    w=multiply(dot(2,pi),f)
# true_amplitude.m:35
    
    K=dot(norm_cons,dig_gain)
# true_amplitude.m:37
    if exist('zpk','builtin'):
        G=zpk(zes,poles,K)
# true_amplitude.m:39
        num,den=tfdata(G,'v',nargout=2)
# true_amplitude.m:40
        h=freqs(num,den,w)
# true_amplitude.m:41
    else:
        if exist('ts_response'):
            G=ts_response('PZ','nm','counts',K,poles,zes)
# true_amplitude.m:43
            h=G.response(f)
# true_amplitude.m:44
        else:
            error('You need either the "zpk" function (Signal Processing Toolbox) or TiSKit')
    
    #figure,plot(log10(f),log10(abs(h)))
## Calibration gain
    
    # find closest sample to reference frequency at which
# the gain should be equal to 1/g
    
    # example: Suppose the calibration gain is 1/2 nm/c at reference period 1s
# that means for frequence 1 Hz, |T(w_ref)|=2 c/nm=1/gc
    
    # |Tnew(w_ref)|=c*|Told(w_ref)|=1/gc  =>  c=1/[|Told(w_ref)|*gc]
#
    scratch,f_ech=min(abs(f - f_ref),nargout=2)
# true_amplitude.m:60
    c=1 / (multiply(cali_gain,abs(h(f_ech))))
# true_amplitude.m:61
    h=dot(c,h)
# true_amplitude.m:62
    ## Adding conjugate mirror to frequency response function to allow division
    
    h=concat([h,fliplr(conj(h))])
# true_amplitude.m:67
    h[1]=1
# true_amplitude.m:68
    
    h[end()]=1
# true_amplitude.m:69
    ## Spectrum division or deconvolution
    
    h=h.T
# true_amplitude.m:74
    B=repmat(h,1,n)
# true_amplitude.m:75
    fr=Y / B
# true_amplitude.m:76
    output=real(ifft(fr,NFFT))
# true_amplitude.m:77
    output=output(arange(1,m),arange())
# true_amplitude.m:78
    return output
    
if __name__ == '__main__':
    pass
    