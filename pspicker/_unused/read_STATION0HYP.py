# Generated with SMOP  0.41-beta
from smop.libsmop import *
# read_STATION0HYP.m

    ### Function amde to read STATION0.HYP file
# Input:    station_file > Name of station file
# Output:   STA > Cell with name,lon,lat(decimal),depth(m)
    
    
@function
def read_STATION0HYP(station_file=None,*args,**kwargs):
    varargin = read_STATION0HYP.varargin
    nargin = read_STATION0HYP.nargin

    #station_file='/Volumes/donnees/SE_DATA/DAT/STATION0.HYP';
    
    ### Read station_file
    
    foc=fopen(station_file,'rt')
# read_STATION0HYP.m:11
    k=0
# read_STATION0HYP.m:12
    count=0
# read_STATION0HYP.m:13
    while logical_not(feof(foc)):

        line=fgetl(foc)
# read_STATION0HYP.m:15
        if size(line) < 10:
            count=count + 1
# read_STATION0HYP.m:17
            if count == 2:
                return STA
            else:
                continue
        else:
            k=k + 1
# read_STATION0HYP.m:24
            if strcmp(line(1),'-'):
                signedepth=- 1
# read_STATION0HYP.m:26
            else:
                signedepth=1
# read_STATION0HYP.m:28
            if strcmp(line(14),'S'):
                signelat=- 1
# read_STATION0HYP.m:31
            else:
                signelat=1
# read_STATION0HYP.m:33
            if strcmp(line(23),'W'):
                signelon=- 1
# read_STATION0HYP.m:36
            else:
                signelon=1
# read_STATION0HYP.m:38
            if strcmp(line(11),'.'):
                lat=str2num(line(arange(7,8))) + str2num(line(arange(9,13))) / 60
# read_STATION0HYP.m:42
            else:
                lat=str2num(line(arange(7,8))) + str2num(concat([line(arange(9,10)),'.',line(arange(11,13))])) / 60
# read_STATION0HYP.m:44
            if strcmp(line(20),'.'):
                lon=str2num(line(arange(15,17))) + str2num(line(arange(18,22))) / 60
# read_STATION0HYP.m:47
            else:
                lon=str2num(line(arange(15,17))) + str2num(concat([line(arange(18,19)),'.',line(arange(20,22))])) / 60
# read_STATION0HYP.m:49
            STA[k,1]=cellarray([strtrim(line(arange(2,6)))])
# read_STATION0HYP.m:52
            STA[k,3]=cellarray([dot(signelat,lat)])
# read_STATION0HYP.m:53
            STA[k,2]=cellarray([dot(signelon,lon)])
# read_STATION0HYP.m:54
            depth=str2num(line(arange(24,27))) / 1000
# read_STATION0HYP.m:55
            STA[k,4]=cellarray([dot(signedepth,depth)])
# read_STATION0HYP.m:56

    
    fclose(foc)
    return STA
    
if __name__ == '__main__':
    pass
    