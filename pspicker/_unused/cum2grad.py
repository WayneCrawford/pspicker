# Generated with SMOP  0.41-beta
from smop.libsmop import *
# cum2grad.m

    # Transform cumulative into function that is shift at 0 for all maxima
    
    
@function
def cum2grad(f_in=None,option=None,*args,**kwargs):
    varargin = cum2grad.varargin
    nargin = cum2grad.nargin

    if isempty(f_in):
        error('f_in is empty!')
    
    tycalpha=zeros(size(f_in))
# cum2grad.m:7
    #tycgamma=zeros(length(f_in),1);
    tikx=loca_ext(1,length(f_in),f_in,'maxi')
# cum2grad.m:9
    tikx=tikx(arange(),1)
# cum2grad.m:10
    clear('a','b')
    a=find(isnan(f_in) == 0,1,'first')
# cum2grad.m:12
    b=find(isnan(f_in) == 0,1,'last')
# cum2grad.m:13
    if isempty(tikx):
        tikx=concat([[a],[b]])
# cum2grad.m:15
    else:
        tikx=concat([[a],[tikx],[b]])
# cum2grad.m:17
    
    tiky=f_in(tikx)
# cum2grad.m:20
    for j in arange(length(tikx) - 1,1,- 1).reshape(-1):
        tycalpha[arange(tikx(j),tikx(j + 1))]=tiky(j + 1)
# cum2grad.m:23
    
    f_out=f_in - tycalpha
# cum2grad.m:26
    f_out[f_out > 0]=0
# cum2grad.m:27
    if strcmp(option,'normalize'):
        f_out=f_out / abs(min(f_out))
# cum2grad.m:29
    
    return f_out
    
if __name__ == '__main__':
    pass
    