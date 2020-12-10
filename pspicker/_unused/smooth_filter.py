# Generated with SMOP  0.41-beta
from smop.libsmop import *
# smooth_filter.m

    # This function is made to smooth traces sort by rows. The smooth use the
# filter function with the advantage that if the trace begins with NaN
# values it computes the smoothing without erasing the first non-nan values
# as would do the simple filter function. 
# Instead of Smooth(1)=[f(i)+f(i-1)..+f(i-n)]/(n+1) with f(i-1) and below
# NaN values we have now Smooth(i)=f(i)/1 Smooth(i+1)=[f(i+1)+f(i)]/2 ...      
# (See doc filter for more details)
# Input:    'M_in' matrix or array containing the data, column sorted
#           'n_smooth' scalar,size of the smoothing window       
# Ouptut:   'M_out' smoothed data
#
    
    
@function
def smooth_filter(M_in=None,n_smooth=None,*args,**kwargs):
    varargin = smooth_filter.varargin
    nargin = smooth_filter.nargin

    m=size(M_in,1)
# smooth_filter.m:15
    n=size(M_in,2)
# smooth_filter.m:16
    l=[]
# smooth_filter.m:18
    for j in arange(1,n).reshape(-1):
        a=find(logical_not(isnan(M_in(arange(),j))),1,'first')
# smooth_filter.m:20
        if isempty(a):
            l=concat([l,j])
# smooth_filter.m:22
    
    M_in[arange(),l]=[]
# smooth_filter.m:26
    m=size(M_in,1)
# smooth_filter.m:27
    n=size(M_in,2)
# smooth_filter.m:28
    M_inter=copy(M_in)
# smooth_filter.m:30
    M_inter[isnan(M_inter)]=0
# smooth_filter.m:32
    M_out=filter(ones(1,n_smooth) / n_smooth,1,M_inter)
# smooth_filter.m:33
    correction_mat=ones(m,n)
# smooth_filter.m:35
    
    for i in arange(1,n).reshape(-1):
        clear('ind_correc')
        ind_correc=find(logical_not(isnan(M_in(arange(),i))),1,'first')
# smooth_filter.m:39
        if isempty(ind_correc):
            continue
        correction_mat[arange(ind_correc,ind_correc + n_smooth - 1),i]=n_smooth / (arange(1,n_smooth)).T
# smooth_filter.m:43
    
    M_out=multiply(M_out,correction_mat)
# smooth_filter.m:46
    M_out[isnan(M_in)]=NaN
# smooth_filter.m:47
    return M_out
    
if __name__ == '__main__':
    pass
    
    