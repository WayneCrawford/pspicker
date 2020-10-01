# Generated with SMOP  0.41-beta
from smop.libsmop import *
# snr2weight.m

    
@function
def snr2weight(input_=None,notation_vec=None,*args,**kwargs):
    varargin = snr2weight.varargin
    nargin = snr2weight.nargin

    thres3=notation_vec(1)
# snr2weight.m:4
    thres2=notation_vec(2)
# snr2weight.m:5
    thres1=notation_vec(3)
# snr2weight.m:6
    thres0=notation_vec(4)
# snr2weight.m:7
    m=length(input_)
# snr2weight.m:8
    note=multiply(ones(m,1),4)
# snr2weight.m:10
    note[input_ <= thres3]=4
# snr2weight.m:12
    note[logical_and((input_ >= thres3),(input_ <= thres2))]=3
# snr2weight.m:13
    note[logical_and((input_ >= thres2),(input_ <= thres1))]=2
# snr2weight.m:14
    note[logical_and((input_ >= thres1),(input_ <= thres0))]=1
# snr2weight.m:15
    note[input_ > thres0]=0
# snr2weight.m:16
    return note
    
if __name__ == '__main__':
    pass
    