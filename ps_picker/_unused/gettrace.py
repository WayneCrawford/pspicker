# Generated with SMOP  0.41-beta
from smop.libsmop import *
# gettrace.m

    
@function
def gettrace(trace=None,char=None,pchar=None,schar=None,*args,**kwargs):
    varargin = gettrace.varargin
    nargin = gettrace.nargin

    # GETTRACE returns traces corresponding to selection criteria
    
    # Usage: [datP,datS,dat_noH]=gettrace(trace,char,pchar,schar)
    
    # Inputs:
    #   trace:  MxN array of seismic traces, one trace per column
    #           The trace length is M, the number of traces is N
    #   char:   N-length character vector containing the component
    #           corresponding to each trace: the only possibilities are
    #           "X", "Y", "Z" and "H"
    #   pchar:  character vector of traces to extract for "P-wave"
    #           data
    #   schar:  character vector of traces to extract for "S-wave"
    #           data
    
    # Outputs:
    #   datP:   traces associated with pchar components
    #   datS:   traces associated with schar components
    #   dat_noH:traces NOT associated with component 'H'
    
    
    # Confirm that the "char" string only has XYZH
    if logical_not(isempty(regexp(char,'[^XYZH]'))):
        error('char array contains a character other than "X","Y","Z" or "H"')
    
    
    # Translate alternative provided component codes to X,Y,or H
    pchar=regexprep(pchar,'3','Z')
# gettrace.m:29
    pchar=regexprep(pchar,'N','Y')
# gettrace.m:30
    pchar=regexprep(pchar,'1','Y')
# gettrace.m:31
    pchar=regexprep(pchar,'E','X')
# gettrace.m:32
    pchar=regexprep(pchar,'2','X')
# gettrace.m:33
    schar=regexprep(schar,'N','Y')
# gettrace.m:34
    schar=regexprep(schar,'E','X')
# gettrace.m:35
    schar=regexprep(schar,'1','Y')
# gettrace.m:36
    schar=regexprep(schar,'2','X')
# gettrace.m:37
    
    ind=[]
# gettrace.m:40
    for i in arange(1,length(pchar)).reshape(-1):
        v=regexpi(char,pchar(i))
# gettrace.m:42
        ind=concat([ind,v])
# gettrace.m:43
    
    #disp(ind)
    datP=trace(arange(),ind)
# gettrace.m:46
    
    ind=[]
# gettrace.m:49
    for i in arange(1,length(schar)).reshape(-1):
        v=regexpi(char,schar(i))
# gettrace.m:51
        ind=concat([ind,v])
# gettrace.m:52
    
    datS=trace(arange(),ind)
# gettrace.m:54
    
    ind=regexpi(char,'[^H]')
# gettrace.m:57
    dat_noH=trace(arange(),ind)
# gettrace.m:58
    return datP,datS,dat_noH
    
if __name__ == '__main__':
    pass
    