# Generated with SMOP  0.41-beta
from smop.libsmop import *
# f_segment.m

   
@function
def f_segment(f=None,*args,**kwargs):
    """
    Calculate the segment between first and last value of function.
    
    Input and output have the same size.  
    :param function: (traces has to be ordered by column)
    :returns: the segment vector
    """
    # f_segment.m:12
    m=size(f,1)
    n=size(f,2)
    segment=NaN(m,n)
    for i in arange(1,n).reshape(-1):
        clear('a','b','ya','yb','lin')
        a=find(logical_not(isnan(f(arange(),i))),1,'first')
        b=find(logical_not(isnan(f(arange(),i))),1,'last')
        ya=f(a,i)
        yb=f(b,i)
        lin=linspace(ya,yb,b - a)
        segment[arange(),i]=inc_nan(lin,m,a)
    return segment


def inc_nan(small_v=None,n=None,i=None,*args,**kwargs):
    """
   Include a vector in a bigger NaN vector
   
    :param small_v: the vector itself
    :param input: length of NaN vector (1 dim)
    :param n: first sample to include the vector "i",
    :returns: vector which has the same size as NaN "big_v"
    """
    # inc_nan.m:10
    if size(small_v,1) == 1:
        small_v=small_v.T
    
    if length(small_v) > n:
        return big_v
    else:
        b=length(small_v)
        big_v=cat(1,NaN(i - 1,1),small_v,NaN(n - (b + i) + 1,1))
    return big_v
    
if __name__ == '__main__':
    pass
    