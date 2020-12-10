# Generated with SMOP  0.41-beta
from smop.libsmop import *
# filterbutter.m

    # Function to filter data with a butterworth filter
# order: order of the butterworth filter
# f1 and f2: low and high cutoffs frequencies
# fech: sample rate
# trace: raw data to filter. trace can be a matrix of data with each column
# representing a trace. The filter will perform column by column
    
    
@function
def filterbutter(order=None,f1=None,f2=None,fech=None,trace=None,*args,**kwargs):
    varargin = filterbutter.varargin
    nargin = filterbutter.nargin

    
    #tracefil=NaN(length(trace),1);
    #new=trace(~isnan(trace));
    #[b,a]=butter(order,[f1 f2]/(fech/2));
    # tracefil=filtfilt(b,a,trace);
    b,a=private_butter(order,concat([f1,f2]) / (fech / 2),nargout=2)
# filterbutter.m:14
    
    tracefil=private_filtfilt(b,a,trace)
# filterbutter.m:15
    
    
    return tracefil
    
if __name__ == '__main__':
    pass
    