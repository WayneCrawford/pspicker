# Generated with SMOP  0.41-beta
from smop.libsmop import *
# trace2FWkurto.m

    
@function
def trace2FWkurto(f=None,h=None,FB=None,T=None,n_smooth=None,first=None,last=None,*args,**kwargs):
    varargin = trace2FWkurto.varargin
    nargin = trace2FWkurto.nargin

    # trace2FWkurto calculates cumulative kurtosis over bandwidths and windows
    
    # Usage:
#   [mean_M,C]=trace2FWkurto(f,h,FB,T,n_smooth,first,last)
    
    # The FW in the name stands for Frequency Window because it computes for
# several frequency bandwidth and window sizes
    
    # Input:    'f' is the raw trace (works on array, not on matrix)
#           'h' is the sampling frequency
#           'FB' is an n-by-2 matrix containing n frequency band
#                   specifications
#           'T' is an m-by-1 vector of window lengths (seconds)
#           'first' is the first sample of interest
#           'last' is ,the last sample of interest
# Ouput:    'M' is the (m*n)-by-(length(f)) kurtosis cumulative matrix
#           'sum_trace' is the mean value of the kurtosis for all windows
#           and all frequency bandwidth
# Example:  trace2FWkurto(f,100,[5 10;15 20],[1 2 3],10,1,1500)
    
    ### More parameters
    
    T=ravel(T)
# trace2FWkurto.m:26
    n=size(FB,1)
# trace2FWkurto.m:27
    m=size(T,1)
# trace2FWkurto.m:28
    l=dot(m,n)
# trace2FWkurto.m:29
    fri=size(f,2)
# trace2FWkurto.m:30
    dsample=1 / h
# trace2FWkurto.m:31
    nsample=size(f,1)
# trace2FWkurto.m:32
    if last > nsample:
        last=copy(nsample)
# trace2FWkurto.m:34
    
    if first < 1:
        first=1
# trace2FWkurto.m:37
    
    Nwin=floor(T / dsample)
# trace2FWkurto.m:39
    ### Trace filtering
    
    k=1
# trace2FWkurto.m:43
    M=NaN(nsample,size(f,2),dot(n,m))
# trace2FWkurto.m:44
    Mf=copy(M)
# trace2FWkurto.m:45
    B=[]
# trace2FWkurto.m:46
    C=[]
# trace2FWkurto.m:47
    for i in arange(1,n).reshape(-1):
        A=filterbutter(3,FB(i,1),FB(i,2),h,f)
# trace2FWkurto.m:49
        A=same_inc(A,first,last)
# trace2FWkurto.m:50
        B=concat([B,A])
# trace2FWkurto.m:51
    
    for j in arange(1,m).reshape(-1):
        # WCC: make sure that the kurtosis window isn't longer than the data window
        if Nwin(j) > (last - first):
            warning('Window length (%2gs) > selected data window (%4gs), skipping!',Nwin(j) / h,(last - first) / h)
            #    warning('Window length (#2gs) > 1/2 selected data window (#4gs)',Nwin(j)/h,(last-first)/h);
        else:
            kurtos=fast_kurtosis(B,Nwin(j))
# trace2FWkurto.m:60
            f_kurto_cum=smooth_filter(kurtos,n_smooth)
# trace2FWkurto.m:61
            kurto_cum=f_cumul(f_kurto_cum)
# trace2FWkurto.m:62
            line=f_segment(kurto_cum)
# trace2FWkurto.m:63
            corr_kurto_cum=kurto_cum - line
# trace2FWkurto.m:64
            C=concat([C,corr_kurto_cum])
# trace2FWkurto.m:65
    
    mean_M=sum(C,2) / (dot(l,fri))
# trace2FWkurto.m:68
    return mean_M,C
    
if __name__ == '__main__':
    pass
    