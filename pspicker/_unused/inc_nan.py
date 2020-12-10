# Generated with SMOP  0.41-beta
from smop.libsmop import *
# inc_nan.m

    #####
# function to include a vector in a bigger NaN vector
# input: length of NaN vector (1 dim) "n", first sample to include the vector "i",
# the vector itself "small_v"
# ouput: vector which has the same size as NaN "big_v"
    
    
@function
def inc_nan(small_v=None,n=None,i=None,*args,**kwargs):
    varargin = inc_nan.varargin
    nargin = inc_nan.nargin

    if size(small_v,1) == 1:
        small_v=small_v.T
# inc_nan.m:10
    
    if length(small_v) > n:
        return big_v
    else:
        b=length(small_v)
# inc_nan.m:16
        big_v=cat(1,NaN(i - 1,1),small_v,NaN(n - (b + i) + 1,1))
# inc_nan.m:17
    
    return big_v
    
if __name__ == '__main__':
    pass
    