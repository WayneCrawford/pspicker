# Generated with SMOP  0.41-beta
from smop.libsmop import *
# cells2amp.m

    
@function
def cells2amp(C1=None,C2=None,*args,**kwargs):
    varargin = cells2amp.varargin
    nargin = cells2amp.nargin

    # Return the cell with the largest amplitudes for each station
    # 
    # Input:
    # 	Cell arrays with the (very ugly) structure:
    #       {1}: station names
    #       {2}: phase identifier???  (int)
    #       {3}: SNRs
    #       {4}: [period, amplitude, index]s
    #       {5}: channel names to write result to
    
    # Output:
    #   Cell array with elements {station, index, amplitude, period}
    
    # Concatenate the cells
    for i in arange(1,size(C1,2)).reshape(-1):
        A[i]=cellarray([concat([[C1[i]],[C2[i]]])])
# cells2amp.m:17
    
    k=1
# cells2amp.m:20
    j=1
# cells2amp.m:21
    if isempty(A[1]):
        Cellout=copy(A)
# cells2amp.m:24
    else:
        while logical_not(isempty(A[1])):

            clear('B','d','M','E','F')
            # Find index of cells whose station name matches the first one
            B=strcmp(A[1],A[1](1))
# cells2amp.m:29
            d=find(B == 1)
# cells2amp.m:30
            M=cell2mat(A[end() - 1])
# cells2amp.m:32
            E=M(d,2)
# cells2amp.m:33
            F=M(M(arange(),2) == max(E),arange())
# cells2amp.m:35
            F=F(1,arange())
# cells2amp.m:36
            G[1][j,1]=A[1](1)
# cells2amp.m:38
            G[2][j,1]=F(3)
# cells2amp.m:39
            G[3][j,1]=F(2)
# cells2amp.m:40
            G[4][j,1]=F(1)
# cells2amp.m:41
            G[5][j,1]=A[5](1)
# cells2amp.m:42
            A=rm_cell_line(d,A)
# cells2amp.m:43
            A[1]
            j=j + 1
# cells2amp.m:45

        Cellout=copy(G)
# cells2amp.m:47
    
    return Cellout
    
if __name__ == '__main__':
    pass
    
    
    