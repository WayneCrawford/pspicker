# Generated with SMOP  0.41-beta
from smop.libsmop import *
# read_selec.m

    # Function made to read auto.out,collect.out, select.out or simple S-files
# and make a structure 
# Input: select.out or collect.out or S-file from seisan SELECT
# Events
#     ->date
#     ->depth
#     ->RMS
#     ->Cell->station;seconds;type;weight
    
    
@function
def read_selec(file=None,*args,**kwargs):
    varargin = read_selec.varargin
    nargin = read_selec.nargin

    fic=fopen(file,'r')
# read_selec.m:12
    count=0
# read_selec.m:13
    clear('Event')
    while logical_not(feof(fic)):

        line=fgetl(fic)
# read_selec.m:17
        if strcmp(line(end()),'1'):
            count=count + 1
# read_selec.m:19
            j=0
# read_selec.m:20
            clear('C','D')
            header=copy(line)
# read_selec.m:22
            Year=str2num(header(arange(2,5)))
# read_selec.m:24
            Month=str2num(header(arange(7,8)))
# read_selec.m:25
            Day=str2num(header(arange(9,10)))
# read_selec.m:26
            Hour=str2num(header(arange(12,13)))
# read_selec.m:27
            Minut=str2num(header(arange(14,15)))
# read_selec.m:28
            Sec=str2num(header(arange(17,20)))
# read_selec.m:29
            Time=datenum(concat([Year,Month,Day,Hour,Minut,Sec]))
# read_selec.m:30
            Coord=concat([str2num(header(arange(24,30))),str2num(header(arange(31,38)))])
# read_selec.m:32
            Depth=str2num(header(arange(39,43)))
# read_selec.m:33
            RMS=str2num(header(arange(52,55)))
# read_selec.m:34
            Mag=str2num(header(arange(56,59)))
# read_selec.m:35
        else:
            if strcmp(line(end()),'6'):
                seed_name=textscan(line,'%s')
# read_selec.m:38
                seed_name=seed_name[1](1)
# read_selec.m:39
            else:
                if strcmp(line(80),' '):
                    data=copy(line)
# read_selec.m:42
                    if strcmp(data(2),' '):
                        Coord=[]
# read_selec.m:44
                        Depth=[]
# read_selec.m:44
                        Mag=[]
# read_selec.m:44
                        RMS=[]
# read_selec.m:44
                        Numberpicks=0
# read_selec.m:44
                        D=cellarray([])
# read_selec.m:44
                    else:
                        while logical_not(strcmp(data(2),' ')):

                            j=j + 1
# read_selec.m:47
                            phour=str2num(data(arange(19,20)))
# read_selec.m:48
                            pmin=str2num(data(arange(21,22)))
# read_selec.m:49
                            psec=str2num(data(arange(23,28)))
# read_selec.m:50
                            ptime=datenum(concat([0,0,0,phour,pmin,psec]))
# read_selec.m:51
                            pweight=str2num(data(arange(15,15)))
# read_selec.m:52
                            pstat=strtrim(data(arange(2,6)))
# read_selec.m:53
                            pphase=data(arange(11,11))
# read_selec.m:54
                            C[j,arange()]=cellarray([cellarray([pstat]),cellarray([pphase]),cellarray([pweight]),cellarray([ptime])])
# read_selec.m:56
                            Numberpicks=size(C,1)
# read_selec.m:57
                            data=fgetl(fic)
# read_selec.m:58

                        for i in arange(1,size(C,2)).reshape(-1):
                            D[i]=cellarray([vertcat(C[arange(1,end()),i])])
# read_selec.m:61
                    Event[count]=struct('SEED_file',seed_name,'Time',Time,'Coord',Coord,'Depth',Depth,'Mag',Mag,'RMS',RMS,'Number_of_picks',Numberpicks,'Cell',cellarray([D]))
# read_selec.m:64

    