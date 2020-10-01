# Generated with SMOP  0.41-beta
from smop.libsmop import *
# loca_ext.m

    
@function
def loca_ext(trace, start, end, type_):
    """
    Return first local extrema
    
    :param trace: waveform data trace
    :param start: start_time
    :param end: end_time
    :param type_: 'maxi' or 'mini'
    """
    y4 = diff(sign(diff(trace.data)))
    if first_ind < 1:
        first_ind == 1    
    if ind_right > length(y4):
        ind_right=length(y4)
    if first_ind > ind_right:
        first_ind=copy(ind_right)
    y4=same_inc(y4,first_ind,ind_right)
# loca_ext.m:10
    if strcmp(type_,'maxi'):
        max_loc=concat([[0],[y4 < 0],[0]])
# loca_ext.m:14
        Fmax=find(max_loc == 1)
# loca_ext.m:15
        indice=copy(Fmax)
# loca_ext.m:16
    else:
        min_loc=concat([[0],[y4 > 0],[0]])
# loca_ext.m:18
        Fmin=find(min_loc == 1)
# loca_ext.m:19
        indice=copy(Fmin)
# loca_ext.m:20
    
    value=f(indice)
# loca_ext.m:23
    A=concat([indice,value])
# loca_ext.m:24
    return A
    
if __name__ == '__main__':
    pass
    