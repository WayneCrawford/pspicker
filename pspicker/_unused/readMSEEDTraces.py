# Generated with SMOP  0.41-beta
from smop.libsmop import *
# readMSEEDTraces.m

    
@function
def readMSEEDTraces(filename=None,ip=None,*args,**kwargs):
    """
    fills traces from rdmseed output
    #   
    # Inputs:
    #   X=first output from rdmseed
    #   I=second output from rdmseed
    #   ip=mappings to use for miniSEED channel names 
    #       .compZ: mapping to Z channel
    #       .compE: mapping to E channel
    #       .compN: mapping to N channel
    #       .compH: mapping to H channel
    
    # Outputs:
    #   stations : list of stations (in same order as in TR*)
    #   traces: cell array of traces, each with the following attributes:
    #           .data:            the data
    #           .sr:              sampling rate
    #           .starttime:       start time (Matlab datenum)
    #           .channelFullName: channel full name (net:sta:cha:comp)
    #   [TRZ,TRN,TRE,TRH]: cell arrays of the indices corresponding to that
    #           channel for each station
    #   P_write,S_write: names of channels to write P and wave picks  to,
    #                   for each station
    
    # Note on mappings:
    #   Mappings are regular expressions that must match the end of the component name.
    #   For example:
    #       '[XE2]'         accepts the letter 'X', 'E', or '2'
    #       '(X|E|2)'       "        "      "   "    "       "
    #       'SHX'           accepts the sequence 'SHX'
    #       '(SHX|SHE|SH2)'	accepts the sequence 'SHX', 'SHE' or 'SH2'
    #       ''              will not accept any component names
    """
    verbose=1
# readMSEEDTraces.m:34
    X,I=rdmseed(filename,nargout=2)
# readMSEEDTraces.m:35
    
    # FROM RDMSEED HELP
#   	- to extract station component n from a multiplexed file:
#  		[X,I] = rdmseed(f);
#  		k = I(n).XBlockIndex;
#  		plot(cat(1,X(k).t),cat(1,X(k).d))
#  		datetick('x')
#  		title(I(n).ChannelFullName)
#  
#         channel_names=Info.ChannelFullName)
    
    stations=unique(cellarray([X.StationIdentifierCode]))
# readMSEEDTraces.m:46
    channels=cellarray([I.ChannelFullName])
# readMSEEDTraces.m:47
    chInds=cellarray([I.XBlockIndex])
# readMSEEDTraces.m:48
    
    traces=cell(length(chInds),1)
# readMSEEDTraces.m:51
    for i in arange(1,length(chInds)).reshape(-1):
        iChs=chInds[i]
# readMSEEDTraces.m:53
        traces[i].data = copy(cat(1,X(iChs).d))
# readMSEEDTraces.m:54
        sr=X(iChs(1)).SampleRate
# readMSEEDTraces.m:55
        channelFullName=channels[i]
# readMSEEDTraces.m:56
        if length(iChs) > 1:
            for j in arange(2,length(iChs)).reshape(-1):
                if X(iChs(j)).SampleRate != sr:
                    error('%s: not all traces have same sample rate!',channelFullName)
                else:
                    startPresent=X(iChs(j)).RecordStartTimeMATLAB
# readMSEEDTraces.m:63
                    startPrevious=X(iChs(j - 1)).RecordStartTimeMATLAB
# readMSEEDTraces.m:64
                    secsPrevious=(length(X(iChs(j - 1)).d)) / sr
# readMSEEDTraces.m:65
                    gap=dot((startPresent - startPrevious),86400.0) - secsPrevious
# readMSEEDTraces.m:66
                    if abs(gap) > 0.1 / sr:
                        outp=copy(channelFullName)
# readMSEEDTraces.m:68
                        outp=concat([outp,sprintf('\nstartPresent = %s\n',datestr(startPresent,'yyyy-mm-ddTHH:MM:SS.FFF'))])
# readMSEEDTraces.m:69
                        outp=concat([outp,sprintf('startPrevious = %s\n',datestr(startPrevious,'yyyy-mm-ddTHH:MM:SS.FFF'))])
# readMSEEDTraces.m:70
                        outp=concat([outp,sprintf('startPresent-startPrevious = %g\n',dot(86400.0,(startPresent - startPrevious)))])
# readMSEEDTraces.m:71
                        outp=concat([outp,sprintf('secsPrevious = %.1f\n',secsPrevious)])
# readMSEEDTraces.m:72
                        outp=concat([outp,sprintf('sr = %g,%d\n',sr)])
# readMSEEDTraces.m:73
                        outp=concat([outp,sprintf('sampsPrev=%d\n',length(X(iChs(j - 1))))])
# readMSEEDTraces.m:74
                        if gap > 0:
                            outp=concat([outp,sprintf('trace offset of %g seconds (%g samples)!',gap,dot(gap,sr))])
# readMSEEDTraces.m:76
                        else:
                            outp=concat([outp,sprintf('trace overlap of %g seconds (%g samples)!',- gap,dot(- gap,sr))])
# readMSEEDTraces.m:78
                        error(outp)
        traces[i].sr = copy(sr)
# readMSEEDTraces.m:85
        traces[i].starttime = copy(X(iChs(1)).RecordStartTimeMATLAB)
# readMSEEDTraces.m:86
        a=strsplit(channelFullName,':')
# readMSEEDTraces.m:87
        traces[i].channelFullName = copy(channelFullName)
# readMSEEDTraces.m:88
        traces[i].network = copy(a[1])
# readMSEEDTraces.m:89
        traces[i].station = copy(a[2])
# readMSEEDTraces.m:90
        traces[i].location = copy(a[3])
# readMSEEDTraces.m:91
        traces[i].channel = copy(a[4])
# readMSEEDTraces.m:92
    
    ## SET UP CELL ARRAYS HOLDING THE INDEXES FOR EACH STATION
    len_sta=length(stations)
# readMSEEDTraces.m:96
    
    TRZ=cell(len_sta,1)
# readMSEEDTraces.m:97
    TRN=cell(len_sta,1)
# readMSEEDTraces.m:98
    TRE=cell(len_sta,1)
# readMSEEDTraces.m:99
    TRH=cell(len_sta,1)
# readMSEEDTraces.m:100
    for i in arange(1,length(stations)).reshape(-1):
        sta_code=strtrim(stations[i])
# readMSEEDTraces.m:103
        TRZ[i]=findChannel(channels,sta_code,ip.compZ)
# readMSEEDTraces.m:104
        TRN[i]=findChannel(channels,sta_code,ip.compN)
# readMSEEDTraces.m:105
        TRE[i]=findChannel(channels,sta_code,ip.compE)
# readMSEEDTraces.m:106
        TRH[i]=findChannel(channels,sta_code,ip.compH)
# readMSEEDTraces.m:107
    
    if 'Z' == ip.putPick_P_Comp:
        P_write=shortChNames(traces,TRZ)
# readMSEEDTraces.m:110
    else:
        if 'N' == ip.putPick_P_Comp:
            P_write=shortChNames(traces,TRN)
# readMSEEDTraces.m:111
        else:
            if 'E' == ip.putPick_P_Comp:
                P_write=shortChNames(traces,TRE)
# readMSEEDTraces.m:112
            else:
                if 'H' == ip.putPick_P_Comp:
                    P_write=shortChNames(traces,TRH)
# readMSEEDTraces.m:113
    
    if 'Z' == ip.putPick_S_Comp:
        S_write=shortChNames(traces,TRZ)
# readMSEEDTraces.m:116
    else:
        if 'N' == ip.putPick_S_Comp:
            S_write=shortChNames(traces,TRN)
# readMSEEDTraces.m:117
        else:
            if 'E' == ip.putPick_S_Comp:
                S_write=shortChNames(traces,TRE)
# readMSEEDTraces.m:118
            else:
                if 'H' == ip.putPick_S_Comp:
                    S_write=shortChNames(traces,TRH)
# readMSEEDTraces.m:119
    
    
    if verbose:
        fprintf('%10s | %-10s | %-10s | %-10s | %-10s |\n','Station','Z','N','E','H')
        disp('-----------+------------+------------+------------+------------|')
        for i in arange(1,length(stations)).reshape(-1):
            fprintf('%10s | %-10s | %-10s | %-10s | %-10s |\n',stations[i],chName(traces,TRZ[i]),chName(traces,TRN[i]),chName(traces,TRE[i]),chName(traces,TRH[i]))
        #         for i=1:length(traces)
#             fprintf('#s #s #10d\n', traces{i}.station,traces{i}.channel,length(traces{i}.data))
#         end
    
    
    
    return stations,traces,TRZ,TRN,TRE,TRH,P_write,S_write
    
if __name__ == '__main__':
    pass
    
    
@function
def shortChNames(traces=None,indexes=None,*args,**kwargs):
    varargin = shortChNames.varargin
    nargin = shortChNames.nargin

    # returns the shortened channel names of the referred-to traces
    
    for i in arange(1,length(indexes)).reshape(-1):
        temp=chName(traces,indexes[i])
# readMSEEDTraces.m:143
        names[i]=temp(concat([1,3]))
# readMSEEDTraces.m:144
    
    return names
    return names
    
if __name__ == '__main__':
    pass
    
    
@function
def chName(traces=None,index=None,*args,**kwargs):
    varargin = chName.varargin
    nargin = chName.nargin

    # returns the channel name of the referred to trace, or returns '---'
    # if the index is empty
    
    if isempty(index):
        name='---'
# readMSEEDTraces.m:155
    else:
        name=traces[index].channel
# readMSEEDTraces.m:157
    
    return name
    return name
    
if __name__ == '__main__':
    pass
    
    
@function
def findChannel(channels=None,station=None,compRegEx=None,*args,**kwargs):
    varargin = findChannel.varargin
    nargin = findChannel.nargin

    # Finds the index in cell array "channels" that matches the station
    # name and component code
    
    # returns 0 if none found
    # error if more than one found
    
    if isempty(compRegEx):
        iMatch=[]
# readMSEEDTraces.m:170
    else:
        searchStr=sprintf('(\\w*):%s:(\\w*):%s$',station,compRegEx)
# readMSEEDTraces.m:172
        iMatch=find(logical_not(cellfun(isempty,regexp(channels,searchStr))))
# readMSEEDTraces.m:174
        if length(iMatch) > 1:
            error('More than one match found for station = "%s" channel="%s"',station,compRegEx)
    
    return iMatch
    return iMatch
    
if __name__ == '__main__':
    pass
    