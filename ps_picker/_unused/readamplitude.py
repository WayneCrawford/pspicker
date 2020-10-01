# Generated with SMOP  0.41-beta
from smop.libsmop import *
# readamplitude.m

    ##### Function that read amplitude.txt file
    
    
@function
def readamplitude(filename=None,*args,**kwargs):
    varargin = readamplitude.varargin
    nargin = readamplitude.nargin

    #filename='amplitude2.txt';
    fic=fopen(filename,'rt')
# readamplitude.m:5
    k=0
# readamplitude.m:7
    while logical_not(feof(fic)):

        line=fgetl(fic)
# readamplitude.m:10
        if strcmp(line(arange(1,4)),'% St'):
            line=fgetl(fic)
# readamplitude.m:12
            while logical_and(logical_not(isempty(line)),(line != - 1)):

                k=k + 1
# readamplitude.m:15
                clear('A')
                A=textscan(line,'%s %s')
# readamplitude.m:17
                stat=A[1][1]
# readamplitude.m:18
                ampli_filename[k,1]=A[2]
# readamplitude.m:19
                stat=sprintf('%s        ',stat)
# readamplitude.m:20
                stations[k,1]=cellarray([stat(arange(1,5))])
# readamplitude.m:21
                line=fgetl(fic)
# readamplitude.m:22

        else:
            if strcmp(line(arange(1,4)),'% Na'):
                line=fgetl(fic)
# readamplitude.m:26
                ampfile=textscan(line,'%s')
# readamplitude.m:27

    
    #A=char(stations)
#stations=mat2cell(stations);
    struc.stations = copy(stations)
# readamplitude.m:34
    struc.amplitudefile = copy(ampli_filename)
# readamplitude.m:35
    fclose(fic)
    return struc
    
if __name__ == '__main__':
    pass
    