# Generated with SMOP  0.41-beta
from smop.libsmop import *
# f_cumul.m

    # Function that calculate the positive gradient cumulative of f
# if the gradient is positive we take the cumulative, if the gradient is
# negative then, as long as the gradient is negative, the value is the one of
# the last positive gradient. The output has then only positive gradients
#      ___/
#     /
# ___/
# Input:    'f' matrix or array containing the data
#      
# Ouptut:   'g' cumulative output matrix
    
    
@function
def f_cumul(f=None,*args,**kwargs):
    varargin = f_cumul.varargin
    nargin = f_cumul.nargin

    m=size(f,1)
# f_cumul.m:14
    n=size(f,2)
# f_cumul.m:15
    inp=copy(f)
# f_cumul.m:16
    inp[isnan(f)]=0
# f_cumul.m:17
    grad_f=diff(inp,1,1)
# f_cumul.m:19
    grad_f=concat([[zeros(1,n)],[grad_f]])
# f_cumul.m:20
    
    grad_f[grad_f < 0]=0
# f_cumul.m:21
    
    g=cumsum(grad_f,1)
# f_cumul.m:23
    
    p=copy(g)
# f_cumul.m:24
    p[isnan(f)]=NaN
# f_cumul.m:25
    corr_mat=ones(1,n)
# f_cumul.m:27
    g[isnan(f)]=NaN
# f_cumul.m:28
    for i in arange(1,n).reshape(-1):
        clear('ind')
        ind=find(logical_not(isnan(g(arange(),i))),1,'first')
# f_cumul.m:32
        if isempty(ind):
            corr_mat[i]=1
# f_cumul.m:34
        else:
            corr_mat[i]=g(ind,i)
# f_cumul.m:36
    
    g=g - dot(ones(m,1),corr_mat)
# f_cumul.m:40
    return g,p
    
if __name__ == '__main__':
    pass
    