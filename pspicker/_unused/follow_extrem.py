# Generated with SMOP  0.41-beta
from smop.libsmop import *
# follow_extrem.m

    
    
@function
def follow_extrem(f, type_='mini', num=2, smooth_vec=[1, 5, 10], option=None, sense=None):
    """
    Find extremas using several smoothed versions of a function
    
    Steps:
      1) locate the 'n' first extremas of the smoothest function
      2) refine these points step by step through less and less smooth functions
    
    :param f: data of interest (unsmoothed)
    :param type:  'mini' or 'maxi', depending on wheter you want to follow minima
        or maxima
    :param num:    number of extrema to follow
    :param smooth_vec: vector containing the different smoothings to apply
    :param option: 'normalize': normalize (what?)
                   anything else: don't normalize
    :param sense: 'first': selects the 'num' first extrema, starting from the right
                  anything else: select the 'num' biggest extrema
    :returns:  value_ext, extrema of unsmoothed function
               ind_ext, indices of extrema
    """    
    # Parameters
    if isempty(f):
        error('f is empty!')
    
    smooth_vec=ravel(smooth_vec)
# follow_extrem.m:29
    m=length(f)
# follow_extrem.m:30
    n=length(smooth_vec)
# follow_extrem.m:31
    smooth_vec=sort(smooth_vec,'descend')
# follow_extrem.m:32
    
    
    M=zeros(m,n)
# follow_extrem.m:36
    Mg=cell(n,1)
# follow_extrem.m:37
    
    
    for i in arange(1,length(smooth_vec)).reshape(-1):
        clear('v','v_smooth')
        v_smooth=smooth_filter(f,smooth_vec(i))
# follow_extrem.m:41
        v_smooth=cum2grad(v_smooth,option)
# follow_extrem.m:42
        Ext_indices=loca_ext(1,length(v_smooth),v_smooth,type_)
# follow_extrem.m:43
        Mg[i]=cellarray([Ext_indices(arange(),1)])
# follow_extrem.m:44
        M[arange(),i]=v_smooth
# follow_extrem.m:45
    
    
    # Choose up to 'num' extrema in the smoothest array
    
    ind_sm=Mg[1]
# follow_extrem.m:50
    ord_sm=M(Mg[1],1)
# follow_extrem.m:51
    L=concat([ind_sm,abs(ord_sm)])
# follow_extrem.m:52
    if strcmp(sense,'first'):
        A=copy(L)
# follow_extrem.m:54
        A[A(arange(),2) < 0.1,arange()]=[]
# follow_extrem.m:55
        if isempty(A):
            b=1
# follow_extrem.m:57
        else:
            A=sortrows(A,- 1)
# follow_extrem.m:59
            b=A(1,1)
# follow_extrem.m:60
            B=A(arange(2,end()),arange())
# follow_extrem.m:61
            if logical_not(isempty(B)):
                B=sortrows(B,- 2)
# follow_extrem.m:63
                try:
                    c=B(arange(1,num - 1),1)
# follow_extrem.m:65
                finally:
                    pass
                b=concat([[b],[c]])
# follow_extrem.m:69
    else:
        L=sortrows(L,- 2)
# follow_extrem.m:73
        try:
            b=L(arange(1,num),1)
# follow_extrem.m:75
        finally:
            pass
    
    if isempty(b):
        ind_ext=[]
# follow_extrem.m:81
        ext2=[]
# follow_extrem.m:82
        return M,ind_ext,ext2
    
    
    ind_ext2[arange(),1]=b
# follow_extrem.m:86
    ext2[arange(),1]=M(b,1)
# follow_extrem.m:88
    
    
    y=zeros(length(b),n)
# follow_extrem.m:92
    for k in arange(2,n).reshape(-1):
        if logical_not(isempty(Mg[k])):
            for j in arange(1,length(b)).reshape(-1):
                clear('differ','ind_differ')
                differ=abs(b(j) - Mg[k])
# follow_extrem.m:98
                if differ > 40:
                    c[j]=b(j)
# follow_extrem.m:100
                    continue
                __,ind_differ=min(differ,nargout=2)
# follow_extrem.m:103
                c[j]=Mg[k](ind_differ(1))
# follow_extrem.m:104
            b=copy(c)
# follow_extrem.m:106
        else:
            warning('No extrema found for %d-sample smoothing of Kurtosis',smooth_vec(n))
        y[arange(),k]=M(b,k)
# follow_extrem.m:110
    
    
    ind_ext[arange(),1]=b
# follow_extrem.m:113
    value_ext[arange(),1]=M(ind_ext,1)
# follow_extrem.m:114
    return M,ind_ext,ext2
    
if __name__ == '__main__':
    pass
    