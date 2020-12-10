# Generated with SMOP  0.41-beta
from smop.libsmop import *
# compare_cell.m

    
@function
def compare_cell(C1=None,C2=None,*args,**kwargs):
    varargin = compare_cell.varargin
    nargin = compare_cell.nargin

    # return indices of cells C1 and C2 whose first elements are the same
    C1[1]
    C2[1]
    
    ind_C1=[]
# compare_cell.m:6
    ind_C2=[]
# compare_cell.m:7
    if logical_not(isempty(C1[1])) and logical_not(isempty(C2[1])):
        for i in arange(1,length(C1[1])).reshape(-1):
            clear('a','k')
            a=strcmp(C1[1](i),C2[1])
# compare_cell.m:12
            k=find(a == 1,1)
# compare_cell.m:13
            if logical_not(isempty(k)):
                ind_C1=concat([[ind_C1],[i]])
# compare_cell.m:15
                ind_C2=concat([[ind_C2],[k]])
# compare_cell.m:16
            else:
                continue
    
    return ind_C1,ind_C2
    
if __name__ == '__main__':
    pass
    
    