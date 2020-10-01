# Generated with SMOP  0.41-beta
from smop.libsmop import *
# same_inc.m
# function that extract part of an initial vector and return a new vector
# which has the same size than the old one filled with NaNs
    
    
@function
def same_inc2(mat_old=None,left=None,right=None,*args,**kwargs):
    varargin = same_inc2.varargin
    nargin = same_inc2.nargin

    m=size(mat_old,1)
# same_inc.m:6
    n=size(mat_old,2)
# same_inc.m:7
    if isempty(mat_old):
        mat_new=[]
# same_inc.m:10
    
    len_=size(mat_old,1)
# same_inc.m:13
    # Change shape of the vector if it's not in column
    
    #formula
    
    try:
        mat_new=cat(1,NaN(left - 1,n),mat_old(arange(left,right),arange()),NaN(len_ - right,n))
# same_inc.m:19
    finally:
        pass
    
    
    return mat_new
    
if __name__ == '__main__':
    pass
    