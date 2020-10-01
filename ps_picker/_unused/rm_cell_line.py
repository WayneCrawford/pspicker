# Generated with SMOP  0.41-beta
from smop.libsmop import *
# rm_cell_line.m

    # rm_cell_line is made to erase lines from a cell ("old_cell"). The line numbers that
# have to be erased are in the vector "vec".
    
    
@function
def rm_cell_line(vec=None,old_cell=None,*args,**kwargs):
    varargin = rm_cell_line.varargin
    nargin = rm_cell_line.nargin

    # Erase elements of a cell array
    
    # indicies are specified in vec
    
    rows=size(old_cell,2)
# rm_cell_line.m:9
    lines=size(old_cell[1],1)
# rm_cell_line.m:10
    if isempty(vec):
        new_cell=copy(old_cell)
# rm_cell_line.m:13
    else:
        if (length(vec) <= lines) and (max(vec) <= lines):
            for i in arange(1,rows).reshape(-1):
                old_cell[i][vec]=[]
# rm_cell_line.m:17
        else:
            disp('Vector of indices to be removed exceed cell dimensions')
        new_cell=copy(old_cell)
# rm_cell_line.m:22
    
    return new_cell
    
if __name__ == '__main__':
    pass
    