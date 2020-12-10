# Generated with SMOP  0.41-beta
from smop.libsmop import *
# private_tukeywin.m

    
@function
def private_tukeywin(L=None,r=None,*args,**kwargs):
    varargin = private_tukeywin.varargin
    nargin = private_tukeywin.nargin

    #PRIVATE_TUKEYWIN tapered cosine window
    
    # Usage: w=private_tukeywin(L,r)
    
    #   returns an L-point Tukey window in the column vector, w. A Tukey 
    #   window is a rectangular window with the first and last r/2 percent 
    #   of the samples equal to parts of a cosine. See Definitions for the 
    #   equation that defines the Tukey window. r is a real number between 
    #   0 and 1. If you input r = 0, you obtain a rectwin window. If you 
    #   input r = 1, you obtain a hann window. r defaults to 0.5.
    
    if r > 1:
        r=1
# private_tukeywin.m:13
    
    
    w=ones(L,1)
# private_tukeywin.m:15
    if r <= 0:
        return w
    
    
    x=linspace(0,1,L)
# private_tukeywin.m:18
    w[x < (r / 2)]=dot(0.5,(1 + cos(dot(dot(2,pi),(x(x < (r / 2)) - r / 2)) / r)))
# private_tukeywin.m:19
    w[x >= (1 - r / 2)]=dot(0.5,(1 + cos(dot(dot(2,pi),(x(x >= (1 - r / 2)) - 1 + r / 2)) / r)))
# private_tukeywin.m:20
    return w
    
if __name__ == '__main__':
    pass
    