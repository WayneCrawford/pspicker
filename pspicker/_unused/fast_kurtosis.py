# Generated with SMOP  0.41-beta
from smop.libsmop import *
# fast_kurtosis.m

    # function to compute kurtosis really quickly by using "filter" function
# Input:    'inp' matrix or array containing the trace
#           'window_sample' number of sample in the sliding window      
# Ouptut:   'out' Kurtosis of the trace
#
    
    
@function
def fast_kurtosis(inp=None,window_sample=None,*args,**kwargs):
    varargin = fast_kurtosis.varargin
    nargin = fast_kurtosis.nargin

    if window_sample == 1:
        window_sample=2
# fast_kurtosis.m:11
    
    # Taking care of NaNs
    
    f=copy(inp)
# fast_kurtosis.m:16
    inp=inp - mean(inp(logical_not(isnan(inp))),1)
# fast_kurtosis.m:17
    inp[isnan(inp) == 1]=0
# fast_kurtosis.m:18
    # Variables
    
    Nwin=copy(window_sample)
# fast_kurtosis.m:22
    # Compute kurtosis
    
    # A=filter(ones(Nwin,1)/Nwin,1,inp);
# f1=A;
# f2=filter(ones(Nwin,1)/Nwin,1,inp.^2);
# f3=filter(ones(Nwin,1)/Nwin,1,inp.^3);
# f4=filter(ones(Nwin,1)/Nwin,1,inp.^4);
# m_2=f2-A.^2;
# m_4=f4-4.*A.*f3+6.*(A.^2).*f2-4.*f1.*A.^3+A.^4;
    m_2=filter(ones(Nwin,1) / Nwin,1,inp ** 2)
# fast_kurtosis.m:33
    m_4=filter(ones(Nwin,1) / Nwin,1,inp ** 4)
# fast_kurtosis.m:34
    out=m_4 / (m_2 ** 2)
# fast_kurtosis.m:36
    out[isnan(f) == 1]=NaN
# fast_kurtosis.m:38
    for i in arange(1,size(inp,2)).reshape(-1):
        a=find(logical_not(isnan(f(arange(),i))),1,'first')
# fast_kurtosis.m:40
        out[arange(a,a + Nwin - 2),i]=NaN
# fast_kurtosis.m:41
    
    return out,A
    
if __name__ == '__main__':
    pass
    