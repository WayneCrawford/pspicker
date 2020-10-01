# Generated with SMOP  0.41-beta
from smop.libsmop import *
# rdmseed.m

    
@function
def rdmseed(varargin=None,*args,**kwargs):
    varargin = rdmseed.varargin
    nargin = rdmseed.nargin

    #RDMSEED Read miniSEED format file.
#	X = RDMSEED(F) reads file F and returns a M-by-1 structure X containing
#	M blocks ("data records") of a miniSEED file with headers, blockettes, 
#	and data in dedicated fields, in particular, for each data block X(i):
#		         t: time vector (DATENUM format)
#		         d: data vector (double)
#		BLOCKETTES: existing blockettes (substructures)
    
    #	Known blockettes are 100, 500, 1000, 1001 and 2000. Others will be
#	ignored with a warning message.
    
    #	X = RDMSEED(F,ENCODINGFORMAT,WORDORDER,RECORDLENGTH), when file F does 
#	not include the Blockette 1000 (like Seismic Handler outputs), specifies:
#		- ENCODINGFORMAT: FDSN code (see below); default is 10 = Steim-1;
#		- WORDORDER: 1 = big-endian (default), 0 = little-endian;
#		- RECORDLENGTH: must be a power of 2, at least 256 (default is 4096).
#	If the file contains Blockette 1000 (which is mandatory in the SEED 
#	convention...), these 3 arguments are ignored except with 'force' option.
    
    #	X = RDMSEED without input argument opens user interface to select the 
#	file from disk.
    
    #	[X,I] = RDMSEED(...) returns a N-by-1 structure I with N the detected 
#	number of different channels, and the following fields:
#	    ChannelFullName: channel name,
#	        XBlockIndex: channel's vector index into X,
#	         ClockDrift: vector of time interval errors, in seconds,
#	                     between each data block (relative to sampling
#	                     period). This can be compared to "Max Clock Drift"
#	                     value of a Blockette 52.
#	                        = 0 in perfect case
#	                        < 0 tends to overlapping
#	                        > 0 tends to gapping
#	  OverlapBlockIndex: index of blocks (into X) having a significant 
#	                     overlap with previous block (less than 0.5
#	                     sampling period).
#	        OverlapTime: time vector of overlapped blocks (DATENUM format).
#	      GapBlockIndex: index of blocks (into X) having a significant gap
#	                     with next block (more than 0.5 sampling period).
#	            GapTime: time vector of gapped blocks (DATENUM format).
    
    #	RDMSEED(...) without output arguments plots the imported signal by 
#	concatenating all the data records, in one single plot if single channel
#	is detected, or subplots for multiplexed file (limited to 10 channels).
#	Gaps are shown with red stars, overlaps with green circles.
    
    #	[...] = RDMSEED(F,...,'be') forces big-endian reading (overwrites the
#	automatic detection of endianness coding, which fails in some cases).
    
    #	[...] = RDMSEED(F,...,'notc') disable time correction.
    
    #	[...] = RDMSEED(F,...,'plot') forces the plot with output arguments.
    
    #	[...] = RDMSEED(F,...,'v') uses verbose mode (displays additional 
#	information and warnings when necessary). Use 'vv' for extras, 'vvv'
#	for debuging.
    
    #	Some instructions for usage of the returned structure:
#	
#	- to get concatenated time and data vectors from a single-channel file:
#		X = rdmseed(f,'plot');
#		t = cat(1,X.t);
#		d = cat(1,X.d);
    
    #	- to get the list of channels in a multiplexed file:
#		[X,I] = rdmseed(f);
#		char(I.ChannelFullName)
    
    #	- to extract the station component n from a multiplexed file:
#		[X,I] = rdmseed(f);
#		k = I(n).XBlockIndex;
#		plot(cat(1,X(k).t),cat(1,X(k).d))
#		datetick('x')
#		title(I(n).ChannelFullName)
    
    #	Known encoding formats are the following FDSN codes:
#		 0: ASCII
#		 1: 16-bit integer
#		 2: 24-bit integer
#		 3: 32-bit integer
#		 4: IEEE float32
#		 5: IEEE float64
#		10: Steim-1
#		11: Steim-2
#		12: GEOSCOPE 24-bit (untested)
#		13: GEOSCOPE 16/3-bit gain ranged
#		14: GEOSCOPE 16/4-bit gain ranged
#		19: Steim-3 (alpha and untested)
    
    #	See also MKMSEED to export data in miniSEED format.
    
    
    #	Author: Francois Beauducel <beauducel@ipgp.fr>
#		Institut de Physique du Globe de Paris
#	Created: 2010-09-17
#	Updated: 2014-06-29
    
    #	Acknowledgments:
#		Ljupco Jordanovski, Jean-Marie Saurel, Mohamed Boubacar, Jonathan Berger,
#		Shahid Ullah, Wayne Crawford, Constanza Pardo, Sylvie Barbier,
#		Robert Chase, Arnaud Lemarchand.
    
    #	References:
#		IRIS (2010), SEED Reference Manual: SEED Format Version 2.4, May 2010,
#		  IFDSN/IRIS/USGS, http://www.iris.edu
#		Trabant C. (2010), libmseed: the Mini-SEED library, IRIS DMC.
#		Steim J.M. (1994), 'Steim' Compression, Quanterra Inc.
    
    #	History:
#		[2014-06-29]
#			- 24-bit uncompressed format tested (bug correction), thanks to
#			  Arnaud Lemarchand.
#		[2014-05-31]
#			- applies the time correction to StartTime and X.t (if needed).
#			- new option 'notc' to disable time correction.
#			- Geoscope 16/4 format passed real data archive tests.
#			- fixes a problem when plotting multiplexed channels (thanks to 
#			  Robert Chase).
#		[2014-03-14]
#			- Improved endianness automatic detection (see comments).
#			- Accepts mixed little/big endian encoding in a single file.
#			- minor fixes.
#		[2013-10-25]
#			- Due to obsolete syntax of bitcmp(0,N) in R2013b, replaces all
#			  by: 2^N-1 (which is much faster...)
#		[2013-02-15]
#			- Tests also DayOfYear in header to determine automatically 
#			  little-endian coding of the file.
#			- Adds option 'be' to force big-endian reading (overwrites
#			  automatic detection).
#		[2012-12-21]
#			- Adds a verbose mode
#		[2012-04-21]
#			- Correct bug with Steim + little-endian coding
#			  (thanks to Shahid Ullah)
#		[2012-03-21]
#			- Adds IDs for warning messages
#		[2011-11-10]
#			- Correct bug with multiple channel name length (thanks to
#			  Jonathan Berger)
#		[2011-10-27]
#			- Add LocationIdentifier to X.ChannelFullName
#		[2011-10-24]
#			- Validation of IEEE double encoding (with PQL)
#			- Import/plot data even with file integrity problem (like PQL)
#		[2011-07-21]
#			- Validation of ASCII encoding format (logs)
#			- Blockettes are now stored in substructures below a single
#			  field X.BLOCKETTES
#			- Add import of blockettes 500 and 2000
#			- Accept multi-channel files with various data coding
#		[2010-10-16]
#			- Alpha-version of Steim-3 decoding...
#			- Extend output parameters with channel detection
#			- Add gaps and overlaps on plots
#			- Add possibility to force the plot
#		[2010-10-02]
#			- Add the input formats for GEOSCOPE multiplexed old data files
#			- Additional output argument with gap and overlap analysis
#			- Create a plot when no output argument are specified
#			- Optimize script coding (30 times faster STEIM decoding!)
#		[2010-09-28]
#			- Correction of a problem with STEIM-1 nibble 3 decoding (one 
#			  32-bit difference)
#			- Add reading of files without blockette 1000 with additional
#			  input arguments (like Seismic Handler output files).
#			- Uses warning() function instead of fprintf().
    
    #	Copyright (c) 2014, Francois Beauducel, covered by BSD License.
#	All rights reserved.
    
    #	Redistribution and use in source and binary forms, with or without 
#	modification, are permitted provided that the following conditions are 
#	met:
    
    #	   * Redistributions of source code must retain the above copyright 
#	     notice, this list of conditions and the following disclaimer.
#	   * Redistributions in binary form must reproduce the above copyright 
#	     notice, this list of conditions and the following disclaimer in 
#	     the documentation and/or other materials provided with the distribution
#	                           
#	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
#	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#	POSSIBILITY OF SUCH DAMAGE.
    
    if nargin > 6:
        error('Too many input arguments.')
    
    # global variables shared with sub-functions
    global f,fid,offset,le,ef,wo,rl,forcebe,verbose,notc,force
    # default input arguments
    makeplot=0
# rdmseed.m:203
    
    verbose=0
# rdmseed.m:204
    
    forcebe=0
# rdmseed.m:205
    
    ef=10
# rdmseed.m:206
    
    wo=1
# rdmseed.m:207
    
    rl=2 ** 12
# rdmseed.m:208
    
    force=0
# rdmseed.m:209
    
    notc=0
# rdmseed.m:210
    
    if nargin < 1:
        filename,pathname=uigetfile('*','Please select a miniSEED file...',nargout=2)
# rdmseed.m:213
        f=fullfile(pathname,filename)
# rdmseed.m:214
    else:
        f=varargin[1]
# rdmseed.m:216
    
    if logical_not(ischar(f)) or logical_not(exist(f,'file')):
        error('File %s does not exist.',f)
    
    if nargin > 1:
        verbose=any(strcmpi(varargin,'v')) + dot(2,any(strcmpi(varargin,'vv'))) + dot(3,any(strcmpi(varargin,'vvv')))
# rdmseed.m:224
        makeplot=any(strcmpi(varargin,'plot'))
# rdmseed.m:226
        forcebe=any(strcmpi(varargin,'be'))
# rdmseed.m:227
        notc=any(strcmpi(varargin,'notc'))
# rdmseed.m:228
        force=any(strcmpi(varargin,'force'))
# rdmseed.m:229
    
    nargs=(makeplot > 0) + (verbose > 0) + (forcebe > 0) + (notc > 0) + (force > 0)
# rdmseed.m:231
    if nargin > (1 + nargs):
        ef=varargin[2]
# rdmseed.m:235
        if logical_not(isnumeric(ef)) or logical_not(any(ef == concat([arange(0,5),arange(10,19),arange(30,33)]))):
            error('Argument ENCODINGFORMAT must be a valid FDSN code value.')
    
    if nargin > (2 + nargs):
        wo=varargin[3]
# rdmseed.m:242
        if logical_not(isnumeric(wo)) or (wo != 0 and wo != 1):
            error('Argument WORDORDER must be 0 or 1.')
    
    if nargin > (3 + nargs):
        rl=varargin[4]
# rdmseed.m:249
        if logical_not(isnumeric(rl)) or rl < 256 or rem(log(rl) / log(2),1) != 0:
            error('Argument RECORDLENGTH must be a power of 2 and greater or equal to 256.')
    
    if nargout == 0:
        makeplot=1
# rdmseed.m:256
    
    # sensible limits for multiplexed files
    max_channels=20
# rdmseed.m:260
    
    max_channel_label=6
# rdmseed.m:261
    
    # file is opened in Big-Endian encoding (this is encouraged by SEED)
    fid=fopen(f,'rb','ieee-be')
# rdmseed.m:264
    le=0
# rdmseed.m:265
    # --- tests if the header is mini-SEED
# the 7th character must be one of the "data header/quality indicator", usually 'D'
    header=fread(fid,20,'*char')
# rdmseed.m:269
    if logical_not(ismember(header(7),'DRMQ')):
        if ismember(header(7),'VAST'):
            s=' (seems to be a SEED Volume)'
# rdmseed.m:272
        else:
            s=''
# rdmseed.m:274
        error('File is not in mini-SEED format%s. Cannot read it.',s)
    
    i=1
# rdmseed.m:279
    offset=0
# rdmseed.m:280
    while offset >= 0:

        X[i]=read_data_record
# rdmseed.m:283
        i=i + 1
# rdmseed.m:284

    
    fclose(fid)
    if nargout > 0:
        varargout[1]=X
# rdmseed.m:290
    
    # --- analyses data
    if makeplot or nargout > 1:
        # test if the file is multiplexed or a single channel
        un=unique(cellstr(char(X.ChannelFullName)))
# rdmseed.m:297
        nc=numel(un)
# rdmseed.m:298
        for i in arange(1,nc).reshape(-1):
            k=find(strcmp(cellstr(char(X.ChannelFullName)),un[i]))
# rdmseed.m:300
            I(i).ChannelFullName = copy(X(k(1)).ChannelFullName)
# rdmseed.m:301
            I(i).XBlockIndex = copy(k)
# rdmseed.m:302
            I(i).ClockDrift = copy((dot(concat([[diff(cat(1,X(k).RecordStartTimeMATLAB))],[NaN]]),86400) - cat(1,X(k).NumberSamples) / cat(1,X(k).SampleRate)) / cat(1,X(k).NumberSamples))
# rdmseed.m:303
            I(i).OverlapBlockIndex = copy(k(find(multiply(multiply(I(i).ClockDrift,cat(1,X(k).NumberSamples)),cat(1,X(k).SampleRate)) < - 0.5) + 1))
# rdmseed.m:304
            I(i).OverlapTime = copy(cat(1,X(I(i).OverlapBlockIndex).RecordStartTimeMATLAB))
# rdmseed.m:305
            I(i).GapBlockIndex = copy(k(find(multiply(multiply(I(i).ClockDrift,cat(1,X(k).NumberSamples)),cat(1,X(k).SampleRate)) > 0.5) + 1))
# rdmseed.m:306
            I(i).GapTime = copy(cat(1,X(I(i).GapBlockIndex).RecordStartTimeMATLAB))
# rdmseed.m:307
    
    if nargout > 1:
        varargout[2]=I
# rdmseed.m:311
    
    # --- plots the data
    if makeplot:
        figure
        xlim=concat([min(cat(1,X.t)),max(cat(1,X.t))])
# rdmseed.m:319
        rl=unique(cat(1,X.DataRecordSize))
# rdmseed.m:322
        if numel(rl) == 1:
            rl_text=sprintf('%d bytes',rl)
# rdmseed.m:324
        else:
            rl_text=sprintf('%d-%d bytes',min(rl),max(rl))
# rdmseed.m:326
        # test if all data records have the same sampling rate
        sr=unique(cat(1,X.SampleRate))
# rdmseed.m:330
        if numel(sr) == 1:
            sr_text=sprintf('%g Hz',sr)
# rdmseed.m:332
        else:
            sr_text=sprintf('%d # samp. rates',numel(sr))
# rdmseed.m:334
        # test if all data records have the same encoding format
        ef=unique(cellstr(cat(1,X.EncodingFormatName)))
# rdmseed.m:338
        if numel(ef) == 1:
            ef_text=sprintf('%s',ef[arange()])
# rdmseed.m:340
        else:
            ef_text=sprintf('%d different encod. formats',numel(ef))
# rdmseed.m:342
        if nc == 1:
            plot(cat(1,X.t),cat(1,X.d))
            hold('on')
            for i in arange(1,length(I.GapBlockIndex)).reshape(-1):
                plot(I.GapTime(i),X(I.GapBlockIndex(i)).d(1),'*r')
            for i in arange(1,length(I.OverlapBlockIndex)).reshape(-1):
                plot(I.OverlapTime(i),X(I.OverlapBlockIndex(i)).d(1),'og')
            hold('off')
            set(gca,'XLim',xlim)
            datetick('x','keeplimits')
            grid('on')
            xlabel(sprintf('Time\n(%s to %s)',datestr(xlim(1)),datestr(xlim(2))))
            ylabel('Counts')
            title(sprintf('mini-SEED file "%s"\n%s (%d rec. @ %s - %g samp. @ %s - %s)',f,un[1],length(X),rl_text,numel(cat(1,X.d)),sr_text,ef_text),'Interpreter','none')
        else:
            # plot is done only for real data channels...
            if nc > max_channels:
                warning('Plot has been limited to %d channels (over %d). See help to manage multiplexed file.',max_channels,nc)
                nc=copy(max_channels)
# rdmseed.m:367
            for i in arange(1,nc).reshape(-1):
                subplot(dot(nc,2),1,dot(i,2) + (arange(- 1,0)))
                k=I(i).XBlockIndex
# rdmseed.m:371
                if logical_not(any(strcmp('ASCII',cellstr(cat(1,X(k).EncodingFormatName))))):
                    plot(cat(1,X(k).t),cat(1,X(k).d))
                    hold('on')
                    for ii in arange(1,length(I(i).GapBlockIndex)).reshape(-1):
                        plot(I(i).GapTime(ii),X(I(i).GapBlockIndex(ii)).d(1),'*r')
                    for ii in arange(1,length(I(i).OverlapBlockIndex)).reshape(-1):
                        plot(I(i).OverlapTime(ii),X(I(i).OverlapBlockIndex(ii)).d(1),'og')
                    hold('off')
                set(gca,'XLim',xlim,'FontSize',8)
                h=ylabel(un[i],'Interpreter','none')
# rdmseed.m:384
                if nc > max_channel_label:
                    set(gca,'YTick',[])
                    set(h,'Rotation',0,'HorizontalAlignment','right','FontSize',8)
                datetick('x','keeplimits')
                set(gca,'XTickLabel',[])
                grid('on')
                if i == 1:
                    title(sprintf('mini-SEED file "%s"\n%d channels (%d rec. @ %s - %g data - %s - %s)',f,length(un),length(X),rl_text,numel(cat(1,X(k).d)),sr_text,ef_text),'Interpreter','none')
                if i == nc:
                    datetick('x','keeplimits')
                    xlabel(sprintf('Time\n(%s to %s)',datestr(xlim(1)),datestr(xlim(2))))
            v=copy(version)
# rdmseed.m:401
            if str2double(v(1)) >= 7:
                linkaxes(findobj(gcf,'type','axes'),'x')
    
    ##########################################################################
    
@function
def read_data_record(*args,**kwargs):
    varargin = read_data_record.varargin
    nargin = read_data_record.nargin

    # read_data_record uses global variables f, fid, offset, le, ef, wo, rl, 
#	and verbose. It reads a data record and returns a structure D.
    
    global f,fid,offset,le,ef,wo,rl,verbose,notc,force
    fseek(fid,offset,'bof')
    # --- read fixed section of Data Header (48 bytes)
    D.SequenceNumber = copy(fread(fid,6,'*char').T)
# rdmseed.m:419
    D.DataQualityIndicator = copy(fread(fid,1,'*char'))
# rdmseed.m:420
    D.ReservedByte = copy(fread(fid,1,'*char'))
# rdmseed.m:421
    D.StationIdentifierCode = copy(fread(fid,5,'*char').T)
# rdmseed.m:422
    D.LocationIdentifier = copy(fread(fid,2,'*char').T)
# rdmseed.m:423
    D.ChannelIdentifier = copy(fread(fid,3,'*char').T)
# rdmseed.m:424
    D.NetworkCode = copy(fread(fid,2,'*char').T)
# rdmseed.m:425
    D.ChannelFullName = copy(sprintf('%s:%s:%s:%s',deblank(D.NetworkCode),deblank(D.StationIdentifierCode),deblank(D.LocationIdentifier),deblank(D.ChannelIdentifier)))
# rdmseed.m:426
    # Start Time decoding
    D.RecordStartTime,swapflag=readbtime
# rdmseed.m:431
    if swapflag:
        if le:
            machinefmt='ieee-be'
# rdmseed.m:434
            le=0
# rdmseed.m:435
        else:
            machinefmt='ieee-le'
# rdmseed.m:437
            le=1
# rdmseed.m:438
        position=ftell(fid)
# rdmseed.m:440
        fclose(fid)
        fid=fopen(f,'rb',machinefmt)
# rdmseed.m:442
        fseek(fid,position,'bof')
        if verbose > 0:
            warning('RDMSEED:DataIntegrity','Sequence # %s: need to switch file encoding to %s...\n',D.SequenceNumber,machinefmt)
    
    D.NumberSamples = copy(fread(fid,1,'uint16'))
# rdmseed.m:451
    # Sample Rate decoding
    SampleRateFactor=fread(fid,1,'int16')
# rdmseed.m:454
    SampleRateMultiplier=fread(fid,1,'int16')
# rdmseed.m:455
    if SampleRateFactor > 0:
        if SampleRateMultiplier >= 0:
            D.SampleRate = copy(dot(SampleRateFactor,SampleRateMultiplier))
# rdmseed.m:458
        else:
            D.SampleRate = copy(dot(- 1,SampleRateFactor) / SampleRateMultiplier)
# rdmseed.m:460
    else:
        if SampleRateMultiplier >= 0:
            D.SampleRate = copy(dot(- 1,SampleRateMultiplier) / SampleRateFactor)
# rdmseed.m:464
        else:
            D.SampleRate = copy(1 / (dot(SampleRateFactor,SampleRateMultiplier)))
# rdmseed.m:466
    
    D.ActivityFlags = copy(fread(fid,1,'uint8'))
# rdmseed.m:470
    D.IOFlags = copy(fread(fid,1,'uint8'))
# rdmseed.m:471
    D.DataQualityFlags = copy(fread(fid,1,'uint8'))
# rdmseed.m:472
    D.NumberBlockettesFollow = copy(fread(fid,1,'uint8'))
# rdmseed.m:473
    D.TimeCorrection = copy(fread(fid,1,'int32'))
# rdmseed.m:474
    
    D.OffsetBeginData = copy(fread(fid,1,'uint16'))
# rdmseed.m:475
    D.OffsetFirstBlockette = copy(fread(fid,1,'uint16'))
# rdmseed.m:476
    # --- read the blockettes
    OffsetNextBlockette=D.OffsetFirstBlockette
# rdmseed.m:479
    D.BLOCKETTES = copy([])
# rdmseed.m:481
    b2000=0
# rdmseed.m:482
    
    for i in arange(1,D.NumberBlockettesFollow).reshape(-1):
        fseek(fid,offset + OffsetNextBlockette,'bof')
        BlocketteType=fread(fid,1,'uint16')
# rdmseed.m:486
        if 1000 == BlocketteType:
            # BLOCKETTE 1000 = Data Only SEED (8 bytes)
            OffsetNextBlockette=fread(fid,1,'uint16')
# rdmseed.m:492
            D.BLOCKETTES.B1000.EncodingFormat = copy(fread(fid,1,'uint8'))
# rdmseed.m:493
            D.BLOCKETTES.B1000.WordOrder = copy(fread(fid,1,'uint8'))
# rdmseed.m:494
            D.BLOCKETTES.B1000.DataRecordLength = copy(fread(fid,1,'uint8'))
# rdmseed.m:495
            D.BLOCKETTES.B1000.Reserved = copy(fread(fid,1,'uint8'))
# rdmseed.m:496
        else:
            if 1001 == BlocketteType:
                # BLOCKETTE 1001 = Data Extension (8 bytes)
                OffsetNextBlockette=fread(fid,1,'uint16')
# rdmseed.m:500
                D.BLOCKETTES.B1001.TimingQuality = copy(fread(fid,1,'uint8'))
# rdmseed.m:501
                D.BLOCKETTES.B1001.Micro_sec = copy(fread(fid,1,'int8'))
# rdmseed.m:502
                D.BLOCKETTES.B1001.Reserved = copy(fread(fid,1,'uint8'))
# rdmseed.m:503
                D.BLOCKETTES.B1001.FrameCount = copy(fread(fid,1,'uint8'))
# rdmseed.m:504
            else:
                if 100 == BlocketteType:
                    # BLOCKETTE 100 = Sample Rate (12 bytes)
                    OffsetNextBlockette=fread(fid,1,'uint16')
# rdmseed.m:508
                    D.BLOCKETTES.B100.ActualSampleRate = copy(fread(fid,1,'float32'))
# rdmseed.m:509
                    D.BLOCKETTES.B100.Flags = copy(fread(fid,1,'uint8'))
# rdmseed.m:510
                    D.BLOCKETTES.B100.Reserved = copy(fread(fid,1,'uint8'))
# rdmseed.m:511
                else:
                    if 500 == BlocketteType:
                        # BLOCKETTE 500 = Timing (200 bytes)
                        OffsetNextBlockette=fread(fid,1,'uint16')
# rdmseed.m:515
                        D.BLOCKETTES.B500.VCOCorrection = copy(fread(fid,1,'float32'))
# rdmseed.m:516
                        D.BLOCKETTES.B500.TimeOfException = copy(readbtime)
# rdmseed.m:517
                        D.BLOCKETTES.B500.MicroSec = copy(fread(fid,1,'int8'))
# rdmseed.m:518
                        D.BLOCKETTES.B500.ReceptionQuality = copy(fread(fid,1,'uint8'))
# rdmseed.m:519
                        D.BLOCKETTES.B500.ExceptionCount = copy(fread(fid,1,'uint16'))
# rdmseed.m:520
                        D.BLOCKETTES.B500.ExceptionType = copy(fread(fid,16,'*char').T)
# rdmseed.m:521
                        D.BLOCKETTES.B500.ClockModel = copy(fread(fid,32,'*char').T)
# rdmseed.m:522
                        D.BLOCKETTES.B500.ClockStatus = copy(fread(fid,128,'*char').T)
# rdmseed.m:523
                    else:
                        if 2000 == BlocketteType:
                            # BLOCKETTE 2000 = Opaque Data (variable length)
                            b2000=b2000 + 1
# rdmseed.m:527
                            OffsetNextBlockette=fread(fid,1,'uint16')
# rdmseed.m:528
                            BlocketteLength=fread(fid,1,'uint16')
# rdmseed.m:529
                            OffsetOpaqueData=fread(fid,1,'uint16')
# rdmseed.m:530
                            D.BLOCKETTES.B2000(b2000).RecordNumber = copy(fread(fid,1,'uint32'))
# rdmseed.m:531
                            D.BLOCKETTES.B2000(b2000).DataWordOrder = copy(fread(fid,1,'uint8'))
# rdmseed.m:532
                            D.BLOCKETTES.B2000(b2000).Flags = copy(fread(fid,1,'uint8'))
# rdmseed.m:533
                            NumberHeaderFields=fread(fid,1,'uint8')
# rdmseed.m:534
                            HeaderFields=splitfield(fread(fid,OffsetOpaqueData - 15,'*char').T,'~')
# rdmseed.m:535
                            D.BLOCKETTES.B2000(b2000).HeaderFields = copy(HeaderFields(arange(1,NumberHeaderFields)))
# rdmseed.m:536
                            # decoded using appropriate format (e.g., Quanterra Q330)
                            D.BLOCKETTES.B2000(b2000).OpaqueData = copy(fread(fid,BlocketteLength - OffsetOpaqueData,'*char').T)
# rdmseed.m:539
                        else:
                            OffsetNextBlockette=fread(fid,1,'uint16')
# rdmseed.m:542
                            if verbose > 0:
                                warning('RDMSEED:UnknownBlockette','Unknown Blockette number %d (%s)!\n',BlocketteType,D.ChannelFullName)
    
    # --- read the data stream
    fseek(fid,offset + D.OffsetBeginData,'bof')
    if logical_not(force) and isfield(D.BLOCKETTES,'B1000'):
        EncodingFormat=D.BLOCKETTES.B1000.EncodingFormat
# rdmseed.m:556
        WordOrder=D.BLOCKETTES.B1000.WordOrder
# rdmseed.m:557
        D.DataRecordSize = copy(2 ** D.BLOCKETTES.B1000.DataRecordLength)
# rdmseed.m:558
    else:
        EncodingFormat=copy(ef)
# rdmseed.m:560
        WordOrder=copy(wo)
# rdmseed.m:561
        D.DataRecordSize = copy(rl)
# rdmseed.m:562
    
    uncoded=0
# rdmseed.m:565
    D.d = copy(NaN)
# rdmseed.m:567
    D.t = copy(NaN)
# rdmseed.m:568
    if 0 == EncodingFormat:
        # --- decoding format: ASCII text
        D.EncodingFormatName = copy(cellarray(['ASCII']))
# rdmseed.m:574
        D.d = copy(fread(fid,D.DataRecordSize - D.OffsetBeginData,'*char').T)
# rdmseed.m:575
    else:
        if 1 == EncodingFormat:
            # --- decoding format: 16-bit integers
            D.EncodingFormatName = copy(cellarray(['INT16']))
# rdmseed.m:579
            dd=fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData) / 2),'*int16')
# rdmseed.m:580
            if xor(logical_not(WordOrder),le):
                dd=swapbytes(dd)
# rdmseed.m:582
            D.d = copy(dd(arange(1,D.NumberSamples)))
# rdmseed.m:584
        else:
            if 2 == EncodingFormat:
                # --- decoding format: 24-bit integers
                D.EncodingFormatName = copy(cellarray(['INT24']))
# rdmseed.m:588
                dd=fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData) / 3),'bit24=>int32')
# rdmseed.m:589
                if xor(logical_not(WordOrder),le):
                    dd=swapbytes(dd)
# rdmseed.m:591
                D.d = copy(dd(arange(1,D.NumberSamples)))
# rdmseed.m:593
            else:
                if 3 == EncodingFormat:
                    # --- decoding format: 32-bit integers
                    D.EncodingFormatName = copy(cellarray(['INT32']))
# rdmseed.m:597
                    dd=fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData) / 4),'*int32')
# rdmseed.m:598
                    if xor(logical_not(WordOrder),le):
                        dd=swapbytes(dd)
# rdmseed.m:600
                    D.d = copy(dd(arange(1,D.NumberSamples)))
# rdmseed.m:602
                else:
                    if 4 == EncodingFormat:
                        # --- decoding format: IEEE floating point
                        D.EncodingFormatName = copy(cellarray(['FLOAT32']))
# rdmseed.m:606
                        dd=fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData) / 4),'*float')
# rdmseed.m:607
                        if xor(logical_not(WordOrder),le):
                            dd=swapbytes(dd)
# rdmseed.m:609
                        D.d = copy(dd(arange(1,D.NumberSamples)))
# rdmseed.m:611
                    else:
                        if 5 == EncodingFormat:
                            # --- decoding format: IEEE double precision floating point
                            D.EncodingFormatName = copy(cellarray(['FLOAT64']))
# rdmseed.m:615
                            dd=fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData) / 8),'*double')
# rdmseed.m:616
                            if xor(logical_not(WordOrder),le):
                                dd=swapbytes(dd)
# rdmseed.m:618
                            D.d = copy(dd(arange(1,D.NumberSamples)))
# rdmseed.m:620
                        else:
                            if cellarray([10,11,19]) == EncodingFormat:
                                # --- decoding formats: STEIM-1 and STEIM-2 compression
		#    (c) Joseph M. Steim, Quanterra Inc., 1994
                                steim=find(EncodingFormat == concat([10,11,19]))
# rdmseed.m:625
                                D.EncodingFormatName = copy(cellarray([sprintf('STEIM%d',steim)]))
# rdmseed.m:626
                                # -- by F. Beauducel, October 2010 --
                                #	1. loads all data into a single 16xM uint32 array
		#	2. gets all nibbles from the first row splitted into 2-bit values
		#	3. for each possible nibble value, selects (find) and decodes
		#	   (bitsplit) all the corresponding words, and stores results
		#	   in a 4xN (STEIM1) or 7xN (STEIM2) array previously filled with
		#	   NaN's. For STEIM2 with nibbles 2 or 3, decodes also dnib values
		#	   (first 2-bit of the word)
		#	5. reduces this array with non-NaN values only
		#	6. integrates with cumsum
                                # This method is about 30 times faster than a 'C-like' loops coding...
                                frame32=fread(fid,concat([16,(D.DataRecordSize - D.OffsetBeginData) / 64]),'*uint32')
# rdmseed.m:643
                                if xor(logical_not(WordOrder),le):
                                    frame32=swapbytes(frame32)
# rdmseed.m:645
                                # specific processes for STEIM-3
                                if steim == 3:
                                    # first bit = 1 means second differences
                                    SecondDiff=bitshift(frame32(1,arange()),- 31)
# rdmseed.m:651
                                    squeezed=bitand(bitshift(frame32(1,arange()),- 24),127)
# rdmseed.m:653
                                    k=find(bitget(squeezed,7))
# rdmseed.m:654
                                    if logical_not(isempty(k)):
                                        moredata24=bitand(frame32(1,k),16777215)
# rdmseed.m:656
                                        k=find(squeezed == 80)
# rdmseed.m:657
                                        if logical_not(isempty(k)):
                                            frame32[1,k]=hex2dec('15555555')
# rdmseed.m:659
                                        k=find(squeezed == 96)
# rdmseed.m:661
                                        if logical_not(isempty(k)):
                                            frame32[1,k]=hex2dec('2aaaaaaa')
# rdmseed.m:663
                                        k=find(squeezed == 112)
# rdmseed.m:665
                                        if logical_not(isempty(k)):
                                            frame32[1,k]=hex2dec('3fffffff')
# rdmseed.m:667
                                # nibbles is an array of the same size as frame32...
                                nibbles=bitand(bitshift(repmat(frame32(1,arange()),16,1),repmat(arange(- 30,0,2),size(frame32,2),1).T),3)
# rdmseed.m:673
                                x0=bitsign(frame32(2,1),32)
# rdmseed.m:674
                                xn=bitsign(frame32(3,1),32)
# rdmseed.m:675
                                if 1 == steim:
                                    # STEIM-1: 3 cases following the nibbles
                                    ddd=dot(NaN,ones(4,numel(frame32)))
# rdmseed.m:681
                                    k=find(nibbles == 1)
# rdmseed.m:682
                                    if logical_not(isempty(k)):
                                        ddd[arange(1,4),k]=bitsplit(frame32(k),32,8)
# rdmseed.m:684
                                    k=find(nibbles == 2)
# rdmseed.m:686
                                    if logical_not(isempty(k)):
                                        ddd[arange(1,2),k]=bitsplit(frame32(k),32,16)
# rdmseed.m:688
                                    k=find(nibbles == 3)
# rdmseed.m:690
                                    if logical_not(isempty(k)):
                                        ddd[1,k]=bitsign(frame32(k),32)
# rdmseed.m:692
                                else:
                                    if 2 == steim:
                                        # STEIM-2: 7 cases following the nibbles and dnib
                                        ddd=dot(NaN,ones(7,numel(frame32)))
# rdmseed.m:697
                                        k=find(nibbles == 1)
# rdmseed.m:698
                                        if logical_not(isempty(k)):
                                            ddd[arange(1,4),k]=bitsplit(frame32(k),32,8)
# rdmseed.m:700
                                        k=find(nibbles == 2)
# rdmseed.m:702
                                        if logical_not(isempty(k)):
                                            dnib=bitshift(frame32(k),- 30)
# rdmseed.m:704
                                            kk=k(dnib == 1)
# rdmseed.m:705
                                            if logical_not(isempty(kk)):
                                                ddd[1,kk]=bitsign(frame32(kk),30)
# rdmseed.m:707
                                            kk=k(dnib == 2)
# rdmseed.m:709
                                            if logical_not(isempty(kk)):
                                                ddd[arange(1,2),kk]=bitsplit(frame32(kk),30,15)
# rdmseed.m:711
                                            kk=k(dnib == 3)
# rdmseed.m:713
                                            if logical_not(isempty(kk)):
                                                ddd[arange(1,3),kk]=bitsplit(frame32(kk),30,10)
# rdmseed.m:715
                                        k=find(nibbles == 3)
# rdmseed.m:718
                                        if logical_not(isempty(k)):
                                            dnib=bitshift(frame32(k),- 30)
# rdmseed.m:720
                                            kk=k(dnib == 0)
# rdmseed.m:721
                                            if logical_not(isempty(kk)):
                                                ddd[arange(1,5),kk]=bitsplit(frame32(kk),30,6)
# rdmseed.m:723
                                            kk=k(dnib == 1)
# rdmseed.m:725
                                            if logical_not(isempty(kk)):
                                                ddd[arange(1,6),kk]=bitsplit(frame32(kk),30,5)
# rdmseed.m:727
                                            kk=k(dnib == 2)
# rdmseed.m:729
                                            if logical_not(isempty(kk)):
                                                ddd[arange(1,7),kk]=bitsplit(frame32(kk),28,4)
# rdmseed.m:731
                                    else:
                                        if 3 == steim:
                                            # STEIM-3: 7 cases following the nibbles
                                            ddd=dot(NaN,ones(9,numel(frame32)))
# rdmseed.m:737
                                            k=find(nibbles == 0)
# rdmseed.m:738
                                            if logical_not(isempty(k)):
                                                ddd[arange(1,2),k]=bitsplit(frame32(k),32,16)
# rdmseed.m:740
                                            k=find(nibbles == 1)
# rdmseed.m:742
                                            if logical_not(isempty(k)):
                                                ddd[arange(1,4),k]=bitsplit(frame32(k),32,8)
# rdmseed.m:744
                                            k=find(nibbles == 2)
# rdmseed.m:746
                                            if logical_not(isempty(k)):
                                                dnib2=bitshift(frame32(k(arange(2,end(),2))),- 30)
# rdmseed.m:748
                                                w60=bitand(frame32(k(arange(2,end(),2))),1073741823) + bitshift(bitand(frame32(k(arange(1,end(),2))),1073741823),30)
# rdmseed.m:749
                                                kk=find(dnib2 == 0)
# rdmseed.m:751
                                                if logical_not(isempty(kk)):
                                                    ddd[arange(1,5),k(dot(2,kk))]=bitsplit(w60,60,12)
# rdmseed.m:753
                                                kk=find(dnib2 == 1)
# rdmseed.m:755
                                                if logical_not(isempty(kk)):
                                                    ddd[arange(1,3),k(dot(2,kk))]=bitsplit(w60,60,20)
# rdmseed.m:757
                                            k=find(nibbles == 3)
# rdmseed.m:760
                                            if logical_not(isempty(k)):
                                                dnib=bitshift(frame32(k),- 27)
# rdmseed.m:762
                                                kk=k(dnib == 24)
# rdmseed.m:763
                                                if logical_not(isempty(kk)):
                                                    ddd[arange(1,9),kk]=bitsplit(frame32(kk),27,3)
# rdmseed.m:765
                                                kk=k(dnib == 25)
# rdmseed.m:767
                                                if logical_not(isempty(kk)):
                                                    ddd[1,kk]=bitsign(frame32(kk),27)
# rdmseed.m:769
                                                kk=k(dnib > 27)
# rdmseed.m:771
                                                if logical_not(isempty(kk)):
                                                    ddd[1,kk]=bitsign(frame32(kk),29)
# rdmseed.m:773
                                # Little-endian coding: needs to swap bytes
                                if logical_not(WordOrder):
                                    ddd=flipud(ddd)
# rdmseed.m:780
                                dd=ddd(logical_not(isnan(ddd)))
# rdmseed.m:782
                                # controls the number of samples
                                if numel(dd) != D.NumberSamples:
                                    if verbose > 1:
                                        warning('RDMSEED:DataIntegrity','Problem in %s sequence # %s [%s]: number of samples in header (%d) does not equal data (%d).\n',D.EncodingFormatName[arange()],D.SequenceNumber,D.RecordStartTimeISO,D.NumberSamples,numel(dd))
                                    if numel(dd) < D.NumberSamples:
                                        D.NumberSamples = copy(numel(dd))
# rdmseed.m:790
                                # rebuilds the data vector by integrating the differences
                                D.d = copy(cumsum(concat([[x0],[dd(arange(2,D.NumberSamples))]])))
# rdmseed.m:795
                                if D.d(end()) != xn:
                                    warning('RDMSEED:DataIntegrity','Problem in %s sequence # %s [%s]: data integrity check failed, last_data=%d, Xn=%d.\n',D.EncodingFormatName[arange()],D.SequenceNumber,D.RecordStartTimeISO,D.d(end()),xn)
                                # for debug purpose...
                                if verbose > 2:
                                    D.dd = copy(dd)
# rdmseed.m:804
                                    D.nibbles = copy(nibbles)
# rdmseed.m:805
                                    D.x0 = copy(x0)
# rdmseed.m:806
                                    D.xn = copy(xn)
# rdmseed.m:807
                            else:
                                if 12 == EncodingFormat:
                                    # --- decoding format: GEOSCOPE multiplexed 24-bit integer
                                    D.EncodingFormatName = copy(cellarray(['GEOSCOPE24']))
# rdmseed.m:812
                                    dd=fread(fid,(D.DataRecordSize - D.OffsetBeginData) / 3,'bit24=>double')
# rdmseed.m:813
                                    if xor(logical_not(WordOrder),le):
                                        dd=swapbytes(dd)
# rdmseed.m:815
                                    D.d = copy(dd(arange(1,D.NumberSamples)))
# rdmseed.m:817
                                else:
                                    if cellarray([13,14]) == EncodingFormat:
                                        # --- decoding format: GEOSCOPE multiplexed 16/3 and 16/4 bit gain ranged
		#	(13): 16/3-bit (bit 15 is unused)
		#	(14): 16/4-bit
		#	bits 15-12 = 3 or 4-bit gain exponent (positive) 
		#	bits 11-0 = 12-bit mantissa (positive)
		#	=> data = (mantissa - 2048) / 2^gain
                                        geoscope=7 + dot(8,(EncodingFormat == 14))
# rdmseed.m:826
                                        D.EncodingFormatName = copy(cellarray([sprintf('GEOSCOPE16-%d',EncodingFormat - 10)]))
# rdmseed.m:827
                                        dd=fread(fid,(D.DataRecordSize - D.OffsetBeginData) / 2,'*uint16')
# rdmseed.m:828
                                        if xor(logical_not(WordOrder),le):
                                            dd=swapbytes(dd)
# rdmseed.m:830
                                        dd=(double(bitand(dd,2 ** 12 - 1)) - 2 ** 11) / 2.0 ** double(bitand(bitshift(dd,- 12),geoscope))
# rdmseed.m:832
                                        D.d = copy(dd(arange(1,D.NumberSamples)))
# rdmseed.m:833
                                    else:
                                        if 15 == EncodingFormat:
                                            # --- decoding format: US National Network compression
                                            D.EncodingFormatName = copy(cellarray(['USNN']))
# rdmseed.m:837
                                            uncoded=1
# rdmseed.m:838
                                        else:
                                            if 16 == EncodingFormat:
                                                # --- decoding format: CDSN 16-bit gain ranged
                                                D.EncodingFormatName = copy(cellarray(['CDSN']))
# rdmseed.m:842
                                                uncoded=1
# rdmseed.m:843
                                            else:
                                                if 17 == EncodingFormat:
                                                    # --- decoding format: Graefenberg 16-bit gain ranged
                                                    D.EncodingFormatName = copy(cellarray(['GRAEFENBERG']))
# rdmseed.m:847
                                                    uncoded=1
# rdmseed.m:848
                                                else:
                                                    if 18 == EncodingFormat:
                                                        # --- decoding format: IPG - Strasbourg 16-bit gain ranged
                                                        D.EncodingFormatName = copy(cellarray(['IPGS']))
# rdmseed.m:852
                                                        uncoded=1
# rdmseed.m:853
                                                    else:
                                                        if 30 == EncodingFormat:
                                                            # --- decoding format: SRO format
                                                            D.EncodingFormatName = copy(cellarray(['SRO']))
# rdmseed.m:857
                                                            uncoded=1
# rdmseed.m:858
                                                        else:
                                                            if 31 == EncodingFormat:
                                                                # --- decoding format: HGLP format
                                                                D.EncodingFormatName = copy(cellarray(['HGLP']))
# rdmseed.m:862
                                                                uncoded=1
# rdmseed.m:863
                                                            else:
                                                                if 32 == EncodingFormat:
                                                                    # --- decoding format: DWWSSN gain ranged format
                                                                    D.EncodingFormatName = copy(cellarray(['DWWSSN']))
# rdmseed.m:867
                                                                    uncoded=1
# rdmseed.m:868
                                                                else:
                                                                    if 33 == EncodingFormat:
                                                                        # --- decoding format: RSTN 16-bit gain ranged
                                                                        D.EncodingFormatName = copy(cellarray(['RSTN']))
# rdmseed.m:872
                                                                        uncoded=1
# rdmseed.m:873
                                                                    else:
                                                                        D.EncodingFormatName = copy(cellarray([sprintf('** Unknown (%d) **',EncodingFormat)]))
# rdmseed.m:876
                                                                        uncoded=1
# rdmseed.m:877
    
    if uncoded:
        error('Sorry, the encoding format "%s" is not yet implemented.',D.EncodingFormatName)
    
    # Applies time correction (if needed)
    D.RecordStartTimeMATLAB = copy(datenum(double(concat([D.RecordStartTime(1),0,D.RecordStartTime(arange(2,5))]))) + dot((logical_and(logical_not(notc),bitand(D.ActivityFlags,2)) == 0),D.TimeCorrection) / 10000.0 / 86400)
# rdmseed.m:886
    tv=datevec(D.RecordStartTimeMATLAB)
# rdmseed.m:888
    doy=datenum(tv(arange(1,3))) - datenum(tv(1),1,0)
# rdmseed.m:889
    D.RecordStartTime = copy(concat([tv(1),doy,tv(arange(4,5)),round(dot(tv(6),10000.0)) / 10000.0]))
# rdmseed.m:890
    D.RecordStartTimeISO = copy(sprintf('%4d-%03d %02d:%02d:%07.4f',D.RecordStartTime))
# rdmseed.m:891
    D.t = copy(D.RecordStartTimeMATLAB)
# rdmseed.m:893
    # makes the time vector and applies time correction (if needed)
    if EncodingFormat > 0:
        D.t = copy(D.t + (arange(0,(D.NumberSamples - 1))).T / (dot(D.SampleRate,86400)))
# rdmseed.m:897
    
    offset=ftell(fid)
# rdmseed.m:900
    fread(fid,1,'char')
    
    if feof(fid):
        offset=- 1
# rdmseed.m:903
    
    ##########################################################################
    
@function
def splitfield(s=None,d=None,*args,**kwargs):
    varargin = splitfield.varargin
    nargin = splitfield.nargin

    # splitfield(S) splits string S of D-character separated field names
    C=textscan(s,'%s','Delimiter',d)
# rdmseed.m:911
    c=C[1]
# rdmseed.m:912
    ##########################################################################
    
@function
def readbtime(*args,**kwargs):
    varargin = readbtime.varargin
    nargin = readbtime.nargin

    # readbtime reads BTIME structure from current opened file and returns
#	D = [YEAR,DAY,HOUR,MINUTE,SECONDS]
    
    global fid,forcebe
    Year=fread(fid,1,'*uint16')
# rdmseed.m:922
    DayOfYear=fread(fid,1,'*uint16')
# rdmseed.m:923
    Hours=fread(fid,1,'uint8')
# rdmseed.m:924
    Minutes=fread(fid,1,'uint8')
# rdmseed.m:925
    Seconds=fread(fid,1,'uint8')
# rdmseed.m:926
    fseek(fid,1,0)
    
    Seconds0001=fread(fid,1,'*uint16')
# rdmseed.m:928
    # Automatic detection of little/big-endian encoding
# -- by F. Beauducel, March 2014 --
    
    # If the 2-byte day is >= 512, the file is not opened in the correct
# endianness. If the day is 1 or 256, there is a possible byte-swap and we
# need to check also the year; but we need to consider what is a valid year:
# - years from 1801 to 2047 are OK (swapbytes >= 2312)
# - years from 2048 to 2055 are OK (swapbytes <= 1800)
# - year 2056 is ambiguous (swapbytes = 2056)
# - years from 2057 to 2311 are OK (swapbytes >= 2312)
# - year 1799 is ambiguous (swapbytes = 1799)
# - year 1800 is suspicious (swapbytes = 2055)
    
    # Thus, the only cases for which we are 'sure' there is a byte-swap, are:
# - day >= 512
# - (day == 1 or day == 256) and (year < 1799 or year > 2311)
    
    # Note: in IRIS libmseed, the test is only year>2050 or year<1920.
    if logical_not(forcebe) and (DayOfYear >= 512 or (ismember(DayOfYear,concat([1,256])) and (Year > 2311 or Year < 1799))):
        swapflag=1
# rdmseed.m:949
        Year=swapbytes(Year)
# rdmseed.m:950
        DayOfYear=swapbytes(DayOfYear)
# rdmseed.m:951
        Seconds0001=swapbytes(Seconds0001)
# rdmseed.m:952
    else:
        swapflag=0
# rdmseed.m:954
    
    d=concat([double(Year),double(DayOfYear),Hours,Minutes,Seconds + double(Seconds0001) / 10000.0])
# rdmseed.m:956
    ##########################################################################
    
@function
def bitsplit(x=None,b=None,n=None,*args,**kwargs):
    varargin = bitsplit.varargin
    nargin = bitsplit.nargin

    # bitsplit(X,B,N) splits the B-bit number X into signed N-bit array
#	X must be unsigned integer class
#	N ranges from 1 to B
#	B is a multiple of N
    
    sign=repmat((arange(b,n,- n)).T,1,size(x,1))
# rdmseed.m:966
    x=repmat(x.T,b / n,1)
# rdmseed.m:967
    d=double(bitand(bitshift(x,flipud(sign - b)),2 ** n - 1)) - dot(double(bitget(x,sign)),2 ** n)
# rdmseed.m:968
    ##########################################################################
    
@function
def bitsign(x=None,n=None,*args,**kwargs):
    varargin = bitsign.varargin
    nargin = bitsign.nargin

    # bitsign(X,N) returns signed double value from unsigned N-bit number X.
# This is equivalent to bitsplit(X,N,N), but the formula is simplified so
# it is much more efficient
    
    d=double(bitand(x,2 ** n - 1)) - multiply(double(bitget(x,n)),2 ** n)
# rdmseed.m:979