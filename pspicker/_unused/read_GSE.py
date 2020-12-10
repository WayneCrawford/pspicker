# Generated with SMOP  0.41-beta
from smop.libsmop import *
# read_GSE.m

    ### Function made to read GSE file and return parameters in structure
# Input:    name > Name of response file GSE
# Output:   info > Structure giving all parameters
    
    
@function
def read_GSE(name=None,*args,**kwargs):
    varargin = read_GSE.varargin
    nargin = read_GSE.nargin

    fic=fopen(name)
# read_GSE.m:6
    k=0
# read_GSE.m:8
    while logical_not(feof(fic)):

        line=fgetl(fic)
# read_GSE.m:10
        k=k + 1
# read_GSE.m:11
        if 'CAL' == line(arange(1,3)):
            info.station = copy(line(arange(6,10)))
# read_GSE.m:14
            info.type = copy(strtrim(line(arange(21,29))))
# read_GSE.m:15
            info.sens = copy(str2num(line(arange(30,40))))
# read_GSE.m:16
            info.sampling = copy(str2num(line(arange(47,56))))
# read_GSE.m:17
        else:
            if 'PAZ' == line(arange(1,3)):
                info.normal = copy(str2num(line(arange(12,25))))
# read_GSE.m:19
                num_poles=str2num(line(arange(42,43)))
# read_GSE.m:20
                num_zeros=str2num(line(arange(46,47)))
# read_GSE.m:21
                p=1
# read_GSE.m:22
                z=1
# read_GSE.m:23
                poles=zeros(num_poles,1)
# read_GSE.m:24
                zers=zeros(num_zeros,1)
# read_GSE.m:25
                while p <= num_poles:

                    line=fgetl(fic)
# read_GSE.m:27
                    clear('A')
                    A=textscan(line,'%s %s')
# read_GSE.m:29
                    poles[p]=complex(str2double(A[1][1]),str2double(A[2][1]))
# read_GSE.m:30
                    p=p + 1
# read_GSE.m:31

                info.poles = copy(poles)
# read_GSE.m:33
                while z <= num_zeros:

                    line=fgetl(fic)
# read_GSE.m:35
                    A=textscan(line,'%s %s')
# read_GSE.m:36
                    zers[z]=complex(str2double(A[1][1]),str2double(A[2][1]))
# read_GSE.m:37
                    z=z + 1
# read_GSE.m:38

                info.zeros = copy(zers)
# read_GSE.m:40
            else:
                if 'DIG' == line(arange(1,3)):
                    info.amplif = copy(str2num(line(arange(10,23))))
# read_GSE.m:42
                    info.digitizer = copy(strtrim(line(arange(37,42))))
# read_GSE.m:43
                else:
                    continue

    
    fclose(fic)
    return info
    
if __name__ == '__main__':
    pass
    
    
    
    