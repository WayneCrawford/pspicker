# Generated with SMOP  0.41-beta
from smop.libsmop import *
# rela2abs.m

    # rela : n-array of picks
# t_begin: n-array beginning of the trace expressed in matlab serial number 
# rsample: sampling frequency
    
    # NB
    
    
@function
def rela2abs(rela=None,t_begin=None,rsample=None,*args,**kwargs):
    varargin = rela2abs.varargin
    nargin = rela2abs.nargin

    m=length(rela)
# rela2abs.m:10
    # conversion from number of sample to seconds
    
    coef_datenum=dot(dot(24,60),60)
# rela2abs.m:14
    a=(rela / rsample) / coef_datenum
# rela2abs.m:15
    b=datenum(a)
# rela2abs.m:16
    abso=t_begin + b
# rela2abs.m:17
    return abso
    
if __name__ == '__main__':
    pass
    