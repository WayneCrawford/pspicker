# Generated with SMOP  0.41-beta
from smop.libsmop import *
# pk2pk.m

    
    
@function
def pk2pk(pick=None,windowR=None,windowL=None,data=None,rsample=None,*args,**kwargs):
    varargin = pk2pk.varargin
    nargin = pk2pk.nargin

    # PK2PK returns the maximum waveform amplitude and corresponding period
# Usage: Amp=pk2pk(pick,windowR,windowL,data,rsample)
    
    # Inputs: 
#   pick=index of reference sample (usually a pick)
#   windowR = seconds after pick for right end of window
#   windowL = seconds after pick for left end of window
#   data=waveform
#   rsample=sampling rate (sps)
# Output:
#   Amp= [period, amplitude, index]
    
    m=size(data,1)
# pk2pk.m:15
    n=size(data,2)
# pk2pk.m:16
    M=zeros(n,3)
# pk2pk.m:17
    NwinL=round(dot(windowL,rsample))
# pk2pk.m:19
    NwinR=round(dot(windowR,rsample))
# pk2pk.m:20
    for i in arange(1,n).reshape(-1):
        clear('Maxx','Minn','ind','Diff','Val','Vat','Vil')
        Maxx=loca_ext(pick - NwinL,pick + NwinR,data(arange(),i),'maxi')
# pk2pk.m:24
        Minn=loca_ext(pick - NwinL,pick + NwinR,data(arange(),i),'mini')
# pk2pk.m:25
        ind=min(concat([size(Maxx,1),size(Minn,1)]))
# pk2pk.m:26
        Maxx=Maxx(arange(1,ind),arange())
# pk2pk.m:27
        Minn=Minn(arange(1,ind),arange())
# pk2pk.m:28
        Diff=Maxx - Minn
# pk2pk.m:29
        Val=Diff(abs(Diff(arange(),2)) == max(abs(Diff(arange(),2))),arange())
# pk2pk.m:30
        Val=abs(Val)
# pk2pk.m:31
        Val[arange(),1]=dot(2.0,Val(arange(),1)) / rsample
# pk2pk.m:32
        Val[arange(),2]=Val(arange(),2) / 2
# pk2pk.m:33
        Vat=Maxx(abs(Diff(arange(),2)) == max(abs(Diff(arange(),2))),arange())
# pk2pk.m:34
        Vil=concat([Val,Vat(arange(),1)])
# pk2pk.m:35
        M[i,arange()]=Vil
# pk2pk.m:36
    
    Amp=M(M(arange(),2) == max(M(arange(),2)),arange())
# pk2pk.m:39
    return Amp
    
if __name__ == '__main__':
    pass
    
    
    
    