# Generated with SMOP  0.41-beta
from smop.libsmop import *
# private_butter.m

    
@function
def tsp_butter(n=None,Wn=None,varargin=None,*args,**kwargs):
    varargin = tsp_butter.varargin
    nargin = tsp_butter.nargin

    #TSP_BUTTER Butterworth digital and analog filter design.
    #   [B,A] = TSP_BUTTER(N,Wn) designs an Nth order lowpass digital
    #   Butterworth filter and returns the filter coefficients in length
    #   N+1 vectors B (numerator) and A (denominator). The coefficients
    #   are listed in descending powers of z. The cut-off frequency
    #   Wn must be 0.0 < Wn < 1.0, with 1.0 corresponding to
    #   half the sample rate.
    
    #   If Wn is a two-element vector, Wn = [W1 W2], BUTTER returns an
    #   order 2N bandpass filter with passband  W1 < W < W2.
    #   [B,A] = BUTTER(N,Wn,'high') designs a highpass filter.
    #   [B,A] = BUTTER(N,Wn,'stop') is a bandstop filter if Wn = [W1 W2].
    
    #   When used with three left-hand arguments, as in
    #   [Z,P,K] = BUTTER(...), the zeros and poles are returned in
    #   length N column vectors Z and P, and the gain in scalar K.
    
    #   When used with four left-hand arguments, as in
    #   [A,B,C,D] = BUTTER(...), state-space matrices are returned.
    
    #   BUTTER(N,Wn,'s'), BUTTER(N,Wn,'high','s') and BUTTER(N,Wn,'stop','s')
    #   design analog Butterworth filters.  In this case, Wn can be bigger
    #   than 1.0.
    
    #   See also BUTTORD, BESSELF, CHEBY1, CHEBY2, ELLIP, FREQZ, FILTER.
    
    #   Author(s): J.N. Little, 1-14-87
    #   	   J.N. Little, 1-14-88, revised
    #   	   L. Shure, 4-29-88, revised
    #   	   T. Krauss, 3-24-93, revised
    #   Copyright (c) 1988-98 by The MathWorks, Inc.
    #   #Revision: 1.22 #  #Date: 1997/12/02 18:36:25 #
    
    #   References:
    #     [1] T. W. Parks and C. S. Burrus, Digital Filter Design,
    #         John Wiley & Sons, 1987, chapter 7, section 7.3.3.
    
    #	Taken without permission (I confess!) from the Matlab signal
    #	toolbox by Wayne Crawford, 4/19/99
    
    btype,analog,errStr=local_iirchk(Wn,varargin[arange()],nargout=3)
# private_butter.m:42
    error(errStr)
    
    if n > 500:
        error('Filter order too large.')
    
    
    # step 1: get analog, pre-warped frequencies
    if logical_not(analog):
        fs=2
# private_butter.m:51
        u=dot(dot(2,fs),tan(dot(pi,Wn) / fs))
# private_butter.m:52
    else:
        u=copy(Wn)
# private_butter.m:54
    
    
    Bw=[]
# private_butter.m:57
    
    if btype == 1:
        Wn=copy(u)
# private_butter.m:60
    else:
        if btype == 2:
            Bw=u(2) - u(1)
# private_butter.m:62
            Wn=sqrt(dot(u(1),u(2)))
# private_butter.m:63
        else:
            if btype == 3:
                Wn=copy(u)
# private_butter.m:65
            else:
                if btype == 4:
                    Bw=u(2) - u(1)
# private_butter.m:67
                    Wn=sqrt(dot(u(1),u(2)))
# private_butter.m:68
    
    
    # step 3: Get N-th order Butterworth analog lowpass prototype
    z,p,k=local_buttap(n,nargout=3)
# private_butter.m:72
    
    a,b,c,d=zp2ss(z,p,k,nargout=4)
# private_butter.m:75
    
    if btype == 1:
        a,b,c,d=local_lp2lp(a,b,c,d,Wn,nargout=4)
# private_butter.m:79
    else:
        if btype == 2:
            a,b,c,d=local_lp2bp(a,b,c,d,Wn,Bw,nargout=4)
# private_butter.m:82
        else:
            if btype == 3:
                a,b,c,d=local_lp2hp(a,b,c,d,Wn,nargout=4)
# private_butter.m:85
            else:
                if btype == 4:
                    a,b,c,d=local_lp2bs(a,b,c,d,Wn,Bw,nargout=4)
# private_butter.m:88
    
    
    # step 5: Use Bilinear transformation to find discrete equivalent:
    if logical_not(analog):
        a,b,c,d=local_bilinear(a,b,c,d,fs,nargout=4)
# private_butter.m:93
    
    
    if nargout == 4:
        num=copy(a)
# private_butter.m:97
        den=copy(b)
# private_butter.m:98
        z=copy(c)
# private_butter.m:99
        p=copy(d)
# private_butter.m:100
    else:
        # Transform to zero-pole-gain and polynomial forms:
        if nargout == 3:
            z,p,k=ss2zp(a,b,c,d,1,nargout=3)
# private_butter.m:104
            z=local_buttzeros(btype,n,Wn,Bw,analog)
# private_butter.m:105
            num=copy(z)
# private_butter.m:106
            den=copy(p)
# private_butter.m:107
            z=copy(k)
# private_butter.m:108
        else:
            den=poly(a)
# private_butter.m:110
            num=local_buttnum(btype,n,Wn,Bw,analog,den)
# private_butter.m:111
    
    return num,den,z,p
    
if __name__ == '__main__':
    pass
    
    #---------------------------------
    
@function
def local_buttnum(btype=None,n=None,Wn=None,Bw=None,analog=None,den=None,*args,**kwargs):
    varargin = local_buttnum.varargin
    nargin = local_buttnum.nargin

    # This internal function returns more exact numerator vectors
    # for the num/den case.
    # Wn input is two element band edge vector
    if analog:
        if 1 == btype:
            b=concat([zeros(1,n),n ** (- n)])
# private_butter.m:125
            b=real(dot(b,polyval(den,dot(- 1j,0))) / polyval(b,dot(- 1j,0)))
# private_butter.m:126
        else:
            if 2 == btype:
                b=concat([zeros(1,n),Bw ** n,zeros(1,n)])
# private_butter.m:128
                b=real(dot(b,polyval(den,dot(- 1j,Wn))) / polyval(b,dot(- 1j,Wn)))
# private_butter.m:129
            else:
                if 3 == btype:
                    b=concat([1,zeros(1,n)])
# private_butter.m:131
                    b=real(dot(b,den(1)) / b(1))
# private_butter.m:132
                else:
                    if 4 == btype:
                        r=dot(dot(1j,Wn),((- 1) ** (arange(0,dot(2,n) - 1)).T))
# private_butter.m:134
                        b=poly(r)
# private_butter.m:135
                        b=real(dot(b,polyval(den,dot(- 1j,0))) / polyval(b,dot(- 1j,0)))
# private_butter.m:136
    else:
        Wn=dot(2,atan2(Wn,4))
# private_butter.m:139
        if 1 == btype:
            r=- ones(n,1)
# private_butter.m:142
            w=0
# private_butter.m:143
        else:
            if 2 == btype:
                r=concat([[ones(n,1)],[- ones(n,1)]])
# private_butter.m:145
                w=copy(Wn)
# private_butter.m:146
            else:
                if 3 == btype:
                    r=ones(n,1)
# private_butter.m:148
                    w=copy(pi)
# private_butter.m:149
                else:
                    if 4 == btype:
                        r=exp(dot(dot(1j,Wn),((- 1) ** (arange(0,dot(2,n) - 1)).T)))
# private_butter.m:151
                        w=0
# private_butter.m:152
        b=poly(r)
# private_butter.m:154
        kern=exp(dot(dot(- 1j,w),(arange(0,length(b) - 1))))
# private_butter.m:156
        b=real(dot(b,(dot(kern,ravel(den)))) / (dot(kern,ravel(b))))
# private_butter.m:157
    
    return b
    
if __name__ == '__main__':
    pass
    
    
@function
def local_buttzeros(btype=None,n=None,Wn=None,Bw=None,analog=None,*args,**kwargs):
    varargin = local_buttzeros.varargin
    nargin = local_buttzeros.nargin

    # This internal function returns more exact zeros.
    # Wn input is two element band edge vector
    if analog:
        # for lowpass and bandpass, don't include zeros at +Inf or -Inf
        if 1 == btype:
            z=zeros(0,1)
# private_butter.m:168
        else:
            if 2 == btype:
                z=zeros(n,1)
# private_butter.m:170
            else:
                if 3 == btype:
                    z=zeros(n,1)
# private_butter.m:172
                else:
                    if 4 == btype:
                        z=dot(dot(1j,Wn),((- 1) ** (arange(0,dot(2,n) - 1)).T))
# private_butter.m:174
    else:
        Wn=dot(2,atan2(Wn,4))
# private_butter.m:177
        if 1 == btype:
            z=- ones(n,1)
# private_butter.m:180
        else:
            if 2 == btype:
                z=concat([[ones(n,1)],[- ones(n,1)]])
# private_butter.m:182
            else:
                if 3 == btype:
                    z=ones(n,1)
# private_butter.m:184
                else:
                    if 4 == btype:
                        z=exp(dot(dot(1j,Wn),((- 1) ** (arange(0,dot(2,n) - 1)).T)))
# private_butter.m:186
    
    return z
    
if __name__ == '__main__':
    pass
    
    #---------------------------------
    
@function
def local_iirchk(Wn=None,varargin=None,*args,**kwargs):
    varargin = local_iirchk.varargin
    nargin = local_iirchk.nargin

    #IIRCHK  Parameter checking for BUTTER, CHEBY1, CHEBY2, and ELLIP.
    #   [btype,analog,errStr] = iirchk(Wn,varargin) returns the
    #   filter type btype (1=lowpass, 2=bandpss, 3=highpass, 4=bandstop)
    #   and analog flag analog (0=digital, 1=analog) given the edge
    #   frequency Wn (either a one or two element vector) and the
    #   optional arguments in varargin.  The variable arguments are
    #   either empty, a one element cell, or a two element cell.
    
    #   errStr is empty if no errors are detected; otherwise it contains
    #   the error message.  If errStr is not empty, btype and analog
    #   are invalid.
    
    #   Copyright (c) 1988-98 by The MathWorks, Inc.
    # #Revision: 1.4 #
    
    errStr=''
# private_butter.m:208
    
    analog=0
# private_butter.m:211
    
    btype=1
# private_butter.m:212
    
    
    if length(Wn) == 1:
        btype=1
# private_butter.m:215
    else:
        if length(Wn) == 2:
            btype=2
# private_butter.m:217
        else:
            errStr='Wn must be a one or two element vector.'
# private_butter.m:219
            return btype,analog,errStr
    
    
    if length(varargin) > 2:
        errStr='Too many input arguments.'
# private_butter.m:224
        return btype,analog,errStr
    
    
    # Interpret and strip off trailing 's' or 'z' argument:
    if logical_not(isempty(varargin)):
        if 's' == lower(varargin[end()]):
            analog=1
# private_butter.m:232
            varargin[end()]=[]
# private_butter.m:233
        else:
            if 'z' == lower(varargin[end()]):
                analog=0
# private_butter.m:235
                varargin[end()]=[]
# private_butter.m:236
            else:
                if length(varargin) > 1:
                    errStr='Analog flag must be either 'z' or 's'.'
# private_butter.m:239
                    return btype,analog,errStr
    
    
    # At this point, varargin will either be empty, or contain a single
    # band type flag.
    
    if length(varargin) == 1:
        if 'low' == lower(varargin[1]):
            btype=1
# private_butter.m:251
        else:
            if 'bandpass' == lower(varargin[1]):
                btype=2
# private_butter.m:253
            else:
                if 'high' == lower(varargin[1]):
                    btype=3
# private_butter.m:255
                else:
                    if 'stop' == lower(varargin[1]):
                        btype=4
# private_butter.m:257
                    else:
                        if nargin == 2:
                            errStr=concat(['Option string must be one of 'high', 'stop',',' 'low', 'bandpass', 'z' or 's'.'])
# private_butter.m:260
                        else:
                            errStr=concat(['Filter type must be one of 'high', 'stop',',' 'low', or 'bandpass'.'])
# private_butter.m:263
                        return btype,analog,errStr
        if 1 == btype:
            if length(Wn) != 1:
                errStr='For the 'low' filter option, Wn must have 1 element.'
# private_butter.m:271
                return btype,analog,errStr
        else:
            if 2 == btype:
                if length(Wn) != 2:
                    errStr='For the 'bandpass' filter option, Wn must have 2 elements.'
# private_butter.m:276
                    return btype,analog,errStr
            else:
                if 3 == btype:
                    if length(Wn) != 1:
                        errStr='For the 'high' filter option, Wn must have 1 element.'
# private_butter.m:281
                        return btype,analog,errStr
                else:
                    if 4 == btype:
                        if length(Wn) != 2:
                            errStr='For the 'stop' filter option, Wn must have 2 elements.'
# private_butter.m:286
                            return btype,analog,errStr
    
    return btype,analog,errStr
    
if __name__ == '__main__':
    pass
    
    #---------------------------------------
    
@function
def local_bilinear(z=None,p=None,k=None,fs=None,fp=None,fp1=None,*args,**kwargs):
    varargin = local_bilinear.varargin
    nargin = local_bilinear.nargin

    #LOCAL_BILINEAR Bilinear transformation with optional frequency prewarping.
    #   [Zd,Pd,Kd] = BILINEAR(Z,P,K,Fs) converts the s-domain transfer
    #   function specified by Z, P, and K to a z-transform discrete
    #   equivalent obtained from the bilinear transformation:
    
    #      H(z) = H(s) |
    #                  | s = 2*Fs*(z-1)/(z+1)
    
    #   where column vectors Z and P specify the zeros and poles, scalar
    #   K specifies the gain, and Fs is the sample frequency in Hz.
    #   [NUMd,DENd] = BILINEAR(NUM,DEN,Fs), where NUM and DEN are
    #   row vectors containing numerator and denominator transfer
    #   function coefficients, NUM(s)/DEN(s), in descending powers of
    #   s, transforms to z-transform coefficients NUMd(z)/DENd(z).
    #   [Ad,Bd,Cd,Dd] = BILINEAR(A,B,C,D,Fs) is a state-space version.
    #   Each of the above three forms of BILINEAR accepts an optional
    #   additional input argument that specifies prewarping. For example,
    #   [Zd,Pd,Kd] = BILINEAR(Z,P,K,Fs,Fp) applies prewarping before
    #   the bilinear transformation so that the frequency responses
    #   before and after mapping match exactly at frequency point Fp
    #   (match point Fp is specified in Hz).
    
    #   See also IMPINVAR.
    
    #   Author(s): J.N. Little, 4-28-87
    #   	   J.N. Little, 5-5-87, revised
    #   Copyright (c) 1988-98 by The MathWorks, Inc.
    #   #Revision: 1.12 #  #Date: 1997/12/02 18:36:15 #
    
    #   Gene Franklin, Stanford Univ., motivated the state-space
    #   approach to the bilinear transformation.
    
    mn,nn=size(z,nargout=2)
# private_butter.m:327
    md,nd=size(p,nargout=2)
# private_butter.m:328
    if (nd == 1 and nn < 2) and nargout != 4:
        if mn > md:
            error('Numerator cannot be higher order than denominator.')
        if nargin == 5:
            fp=dot(dot(2,pi),fp)
# private_butter.m:335
            fs=fp / tan(fp / fs / 2)
# private_butter.m:336
        else:
            fs=dot(2,fs)
# private_butter.m:338
        z=z(finite(z))
# private_butter.m:340
        pd=(1 + p / fs) / (1 - p / fs)
# private_butter.m:341
        zd=(1 + z / fs) / (1 - z / fs)
# private_butter.m:342
        kd=(dot(k,prod(fs - z)) / prod(fs - p))
# private_butter.m:344
        zd=concat([[zd],[- ones(length(pd) - length(zd),1)]])
# private_butter.m:345
    else:
        if (md == 1 and mn == 1) or nargout == 4:
            if nargout == 4:
                a=copy(z)
# private_butter.m:349
                b=copy(p)
# private_butter.m:349
                c=copy(k)
# private_butter.m:349
                d=copy(fs)
# private_butter.m:349
                fs=copy(fp)
# private_butter.m:349
                error(abcdchk(a,b,c,d))
                if nargin == 6:
                    fp=copy(fp1)
# private_butter.m:352
                    fp=dot(dot(2,pi),fp)
# private_butter.m:353
                    fs=fp / tan(fp / fs / 2) / 2
# private_butter.m:354
            else:
                if nn > nd:
                    error('Numerator cannot be higher order than denominator.')
                num=copy(z)
# private_butter.m:360
                den=copy(p)
# private_butter.m:360
                if nargin == 4:
                    fp=copy(fs)
# private_butter.m:362
                    fs=copy(k)
# private_butter.m:362
                    fp=dot(dot(2,pi),fp)
# private_butter.m:363
                    fs=fp / tan(fp / fs / 2) / 2
# private_butter.m:364
                else:
                    fs=copy(k)
# private_butter.m:366
                # Put num(s)/den(s) in state-space canonical form.
                a,b,c,d=tf2ss(num,den,nargout=4)
# private_butter.m:369
            # Now do state-space version of bilinear transformation:
            t=1 / fs
# private_butter.m:372
            r=sqrt(t)
# private_butter.m:373
            t1=eye(size(a)) + dot(a,t) / 2
# private_butter.m:374
            t2=eye(size(a)) - dot(a,t) / 2
# private_butter.m:375
            ad=numpy.linalg.solve(t2,t1)
# private_butter.m:376
            bd=dot(t / r,(numpy.linalg.solve(t2,b)))
# private_butter.m:377
            cd=dot(r,c) / t2
# private_butter.m:378
            dd=dot(dot(c / t2,b),t) / 2 + d
# private_butter.m:379
            if nargout == 4:
                zd=copy(ad)
# private_butter.m:381
                pd=copy(bd)
# private_butter.m:381
                kd=copy(cd)
# private_butter.m:381
            else:
                # Convert back to transfer function form:
                p=poly(ad)
# private_butter.m:384
                zd=poly(ad - dot(bd,cd)) + dot((dd - 1),p)
# private_butter.m:385
                pd=copy(p)
# private_butter.m:386
        else:
            error('First two arguments must have the same orientation.')
    
    return zd,pd,kd,dd
    
if __name__ == '__main__':
    pass
    
    #---------------------------------------
    
@function
def local_buttap(n=None,*args,**kwargs):
    varargin = local_buttap.varargin
    nargin = local_buttap.nargin

    #BUTTAP Butterworth analog lowpass filter prototype.
    #   [Z,P,K] = BUTTAP(N) returns the zeros, poles, and gain
    #   for an N-th order normalized prototype Butterworth analog
    #   lowpass filter.  The resulting filter has N poles around
    #   the unit circle in the left half plane, and no zeros.
    
    #   See also BUTTER, CHEB1AP, CHEB2AP, ELLIPAP.
    
    #   Author(s): J.N. Little and J.O. Smith, 1-14-87
    #   	   L. Shure, 1-13-88, revised
    #   Copyright (c) 1988-98 by The MathWorks, Inc.
    #   #Revision: 1.11 #  #Date: 1997/11/26 20:13:06 #
    
    # Poles are on the unit circle in the left-half plane.
    z=[]
# private_butter.m:409
    p=exp(dot(1j,(dot(pi,(arange(1,n - 1,2))) / (dot(2,n)) + pi / 2)))
# private_butter.m:410
    p=concat([[p],[conj(p)]])
# private_butter.m:411
    p=ravel(p)
# private_butter.m:412
    if rem(n,2) == 1:
        p=concat([[p],[- 1]])
# private_butter.m:414
    
    k=real(prod(- p))
# private_butter.m:416
    return z,p,k
    
if __name__ == '__main__':
    pass
    
    #----------------------------------------
    
@function
def local_lp2bp(a=None,b=None,c=None,d=None,wo=None,bw=None,*args,**kwargs):
    varargin = local_lp2bp.varargin
    nargin = local_lp2bp.nargin

    #LP2BP Lowpass to bandpass analog filter transformation.
    #   [NUMT,DENT] = LP2BP(NUM,DEN,Wo,Bw) transforms the lowpass filter
    #   prototype NUM(s)/DEN(s) with unity cutoff frequency to a
    #   bandpass filter with center frequency Wo and bandwidth Bw.
    #   [AT,BT,CT,DT] = LP2BP(A,B,C,D,Wo,Bw) does the same when the
    #   filter is described in state-space form.
    
    #   Author(s): J.N. Little and G.F. Franklin, 8-4-87
    #   Copyright (c) 1988-98 by The MathWorks, Inc.
    #   #Revision: 1.11 #  #Date: 1997/12/02 19:16:59 #
    
    if nargin == 4:
        # handle column vector inputs: convert to rows
        if size(a,2) == 1:
            a=ravel(a).T
# private_butter.m:435
        if size(b,2) == 1:
            b=ravel(b).T
# private_butter.m:438
        # Transform to state-space
        wo=copy(c)
# private_butter.m:441
        bw=copy(d)
# private_butter.m:442
        a,b,c,d=tf2ss(a,b,nargout=4)
# private_butter.m:443
    
    
    error(abcdchk(a,b,c,d))
    __,nb=size(b,nargout=2)
# private_butter.m:447
    mc,ma=size(c,nargout=2)
# private_butter.m:448
    
    q=wo / bw
# private_butter.m:451
    at=dot(wo,concat([[a / q,eye(ma)],[- eye(ma),zeros(ma)]]))
# private_butter.m:452
    bt=dot(wo,concat([[b / q],[zeros(ma,nb)]]))
# private_butter.m:453
    ct=concat([c,zeros(mc,ma)])
# private_butter.m:454
    dt=copy(d)
# private_butter.m:455
    if nargin == 4:
        # Transform back to transfer function
        z,k=tzero(at,bt,ct,dt,nargout=2)
# private_butter.m:459
        num=dot(k,poly(z))
# private_butter.m:460
        den=poly(at)
# private_butter.m:461
        at=copy(num)
# private_butter.m:462
        bt=copy(den)
# private_butter.m:463
    
    return at,bt,ct,dt
    
if __name__ == '__main__':
    pass
    
    #------------------------------------
    
@function
def local_lp2bs(a=None,b=None,c=None,d=None,wo=None,bw=None,*args,**kwargs):
    varargin = local_lp2bs.varargin
    nargin = local_lp2bs.nargin

    #LP2BS Lowpass to bandstop analog filter transformation.
    #   [NUMT,DENT] = LP2BS(NUM,DEN,Wo,Bw) transforms the lowpass filter
    #   prototype NUM(s)/DEN(s) with unity cutoff frequency to a
    #   bandstop filter with center frequency Wo and bandwidth Bw.
    #   [AT,BT,CT,DT] = LP2BS(A,B,C,D,Wo,Bw) does the same when the
    #   filter is described in state-space form.
    
    #   Author(s): J.N. Little and G.F. Franklin, 8-4-87
    #   Copyright (c) 1988-98 by The MathWorks, Inc.
    #   #Revision: 1.11 #  #Date: 1997/12/02 18:37:06 #
    
    if nargin == 4:
        # handle column vector inputs: convert to rows
        if size(a,2) == 1:
            a=ravel(a).T
# private_butter.m:483
        if size(b,2) == 1:
            b=ravel(b).T
# private_butter.m:486
        # Transform to state-space
        wo=copy(c)
# private_butter.m:489
        bw=copy(d)
# private_butter.m:490
        a,b,c,d=tf2ss(a,b,nargout=4)
# private_butter.m:491
    
    
    error(abcdchk(a,b,c,d))
    __,nb=size(b,nargout=2)
# private_butter.m:495
    mc,ma=size(c,nargout=2)
# private_butter.m:496
    
    q=wo / bw
# private_butter.m:499
    at=concat([[dot(wo / q,inv(a)),dot(wo,eye(ma))],[dot(- wo,eye(ma)),zeros(ma)]])
# private_butter.m:500
    bt=- concat([[dot(wo / q,(numpy.linalg.solve(a,b)))],[zeros(ma,nb)]])
# private_butter.m:501
    ct=concat([c / a,zeros(mc,ma)])
# private_butter.m:502
    dt=d - dot(c / a,b)
# private_butter.m:503
    if nargin == 4:
        # Transform back to transfer function
        z,k=tzero(at,bt,ct,dt,nargout=2)
# private_butter.m:507
        num=dot(k,poly(z))
# private_butter.m:508
        den=poly(at)
# private_butter.m:509
        at=copy(num)
# private_butter.m:510
        bt=copy(den)
# private_butter.m:511
    
    return at,bt,ct,dt
    
if __name__ == '__main__':
    pass
    
    #---------------------------------------------
    
@function
def local_lp2hp(a=None,b=None,c=None,d=None,wo=None,*args,**kwargs):
    varargin = local_lp2hp.varargin
    nargin = local_lp2hp.nargin

    #LP2HP Lowpass to highpass analog filter transformation.
    #   [NUMT,DENT] = LP2HP(NUM,DEN,Wo) transforms the lowpass filter
    #   prototype NUM(s)/DEN(s) with unity cutoff frequency to a
    #   highpass filter with cutoff frequency Wo.
    #   [AT,BT,CT,DT] = LP2HP(A,B,C,D,Wo) does the same when the
    #   filter is described in state-space form.
    
    #   Author(s): J.N. Little and G.F. Franklin, 8-4-87
    #   Copyright (c) 1988-98 by The MathWorks, Inc.
    #   #Revision: 1.10 #  #Date: 1997/11/26 20:13:03 #
    
    if nargin == 3:
        # handle column vector inputs: convert to rows
        if size(a,2) == 1:
            a=ravel(a).T
# private_butter.m:531
        if size(b,2) == 1:
            b=ravel(b).T
# private_butter.m:534
        # Transform to state-space
        wo=copy(c)
# private_butter.m:537
        a,b,c,d=tf2ss(a,b,nargout=4)
# private_butter.m:538
    
    
    error(abcdchk(a,b,c,d))
    ma,nb=size(b,nargout=2)
# private_butter.m:542
    mc,ma=size(c,nargout=2)
# private_butter.m:543
    
    at=dot(wo,inv(a))
# private_butter.m:546
    bt=dot(- wo,(numpy.linalg.solve(a,b)))
# private_butter.m:547
    ct=c / a
# private_butter.m:548
    dt=d - dot(c / a,b)
# private_butter.m:549
    if nargin == 3:
        # Transform back to transfer function
        z,k=tzero(at,bt,ct,dt,nargout=2)
# private_butter.m:553
        num=dot(k,poly(z))
# private_butter.m:554
        den=poly(at)
# private_butter.m:555
        at=copy(num)
# private_butter.m:556
        bt=copy(den)
# private_butter.m:557
    
    return at,bt,ct,dt
    
if __name__ == '__main__':
    pass
    
    #------------------------------------
    
@function
def local_lp2lp(a=None,b=None,c=None,d=None,wo=None,*args,**kwargs):
    varargin = local_lp2lp.varargin
    nargin = local_lp2lp.nargin

    #LP2LP Lowpass to lowpass analog filter transformation.
    #   [NUMT,DENT] = LP2LP(NUM,DEN,Wo) transforms the lowpass filter
    #   prototype NUM(s)/DEN(s) with unity cutoff frequency to a
    #   lowpass filter with cutoff frequency Wo.
    #   [AT,BT,CT,DT] = LP2LP(A,B,C,D,Wo) does the same when the
    #   filter is described in state-space form.
    
    #   Author(s): J.N. Little and G.F. Franklin, 8-4-87
    #   Copyright (c) 1988-98 by The MathWorks, Inc.
    #   #Revision: 1.10 #  #Date: 1997/11/26 20:13:40 #
    
    if nargin == 3:
        # handle column vector inputs: convert to rows
        if size(a,2) == 1:
            a=ravel(a).T
# private_butter.m:577
        if size(b,2) == 1:
            b=ravel(b).T
# private_butter.m:580
        # Transform to state-space
        wo=copy(c)
# private_butter.m:583
        a,b,c,d=tf2ss(a,b,nargout=4)
# private_butter.m:584
    
    
    error(abcdchk(a,b,c,d))
    ma,nb=size(b,nargout=2)
# private_butter.m:588
    mc,ma=size(c,nargout=2)
# private_butter.m:589
    
    at=dot(wo,a)
# private_butter.m:592
    bt=dot(wo,b)
# private_butter.m:593
    ct=copy(c)
# private_butter.m:594
    dt=copy(d)
# private_butter.m:595
    if nargin == 3:
        # Transform back to transfer function
        z,k=tzero(at,bt,ct,dt,nargout=2)
# private_butter.m:599
        num=dot(k,poly(z))
# private_butter.m:600
        den=poly(at)
# private_butter.m:601
        at=copy(num)
# private_butter.m:602
        bt=copy(den)
# private_butter.m:603
    
    return at,bt,ct,dt
    
if __name__ == '__main__':
    pass
    