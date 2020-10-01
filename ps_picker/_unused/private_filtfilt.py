# Generated with SMOP  0.41-beta
from smop.libsmop import *
# private_filtfilt.m

    
@function
def tsp_filtfilt(b=None,a=None,x=None,*args,**kwargs):
    varargin = tsp_filtfilt.varargin
    nargin = tsp_filtfilt.nargin

    #FILTFILT Zero-phase forward and reverse digital filtering.
    #	Y = FILTFILT(B, A, X) filters the data in vector X with the
    #	filter described by vectors A and B to create the filtered
    #	data Y.  The filter is described by the difference equation:
    
    #	  y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
    #	                   - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
    
    #	After filtering in the forward direction, the filtered
    #	sequence is then reversed and run back through the filter.
    #	The resulting sequence has precisely zero-phase distortion
    #	and double the filter order.  Care is taken to minimize
    #	startup and ending transients by matching initial conditions.
    
    #	The length of the input x must be more than three times
    #	the filter order, defined as max(length(b)-1,length(a)-1).
    
    
    narginchk(3,3)
    if (logical_or(logical_or(isempty(b),isempty(a)),isempty(x))):
        y=[]
# private_filtfilt.m:22
        return y
    
    
    m,n=size(x,nargout=2)
# private_filtfilt.m:26
    #     if (n>1)&(m>1)
#         error('Only works for vector input.')
#     end
    if m == 1:
        x=ravel(x)
# private_filtfilt.m:31
        isRow=1
# private_filtfilt.m:32
        m=copy(n)
# private_filtfilt.m:33
        n=1
# private_filtfilt.m:33
    else:
        isRow=0
# private_filtfilt.m:35
    
    len_=size(x,1)
# private_filtfilt.m:37
    
    b=ravel(b).T
# private_filtfilt.m:38
    a=ravel(a).T
# private_filtfilt.m:39
    nb=length(b)
# private_filtfilt.m:40
    na=length(a)
# private_filtfilt.m:41
    nfilt=max(nb,na)
# private_filtfilt.m:42
    nfact=dot(3,(nfilt - 1))
# private_filtfilt.m:44
    
    
    if (len_ <= nfact):
        error('Data must have length more than 3 times filter order.')
    
    
    # set up initial conditions to remove dc offset problems at the beginning and
    # end of the sequence
    __,zi=filter(b,a,ones(nfact,1),nargout=2)
# private_filtfilt.m:52
    
    # method".  Slopes of original and extrapolated sequences match at
    # the end points.
    # This reduces end effects.
    head=dot(2,x(1)) - x(arange((nfact + 1),2,- 1))
# private_filtfilt.m:58
    head=ravel(head)
# private_filtfilt.m:58
    tail=dot(2,x(len_)) - x(arange((len_ - 1),len_ - nfact,- 1))
# private_filtfilt.m:59
    tail=ravel(tail)
# private_filtfilt.m:59
    y=concat([[dot(head,ones(1,n))],[x],[dot(tail,ones(1,n))]])
# private_filtfilt.m:60
    
    y=filter(b,a,y,concat([dot(zi,y(1))]))
# private_filtfilt.m:63
    y=y(arange(length(y),1,- 1),arange())
# private_filtfilt.m:64
    y=filter(b,a,y,concat([dot(zi,y(1))]))
# private_filtfilt.m:65
    y=y(arange(length(y),1,- 1),arange())
# private_filtfilt.m:66
    
    y[concat([arange(1,nfact),len_ + nfact + (arange(1,nfact))]),arange()]=[]
# private_filtfilt.m:69
    if isRow:
        y=y.T
# private_filtfilt.m:72
    
    return y
    
if __name__ == '__main__':
    pass
    