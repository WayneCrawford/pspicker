"""
Pole-Zero representation of an instrument

Follows SEED Pole-Zero protocol for Analog Stages, radian format
(Seed Manual v2.4, Appendix C)

G(f) = S_d * A_0 * prod(s - z_n) / prod(s - p_n) = S_d * A_0 * Hp(s)
where s=j*2*pi*f

The corresponding attributes are:
    S_d: gain
    A_0: norm_factor
    G(f0): passband_gain, with f0: ref_freq

The user can specify gain or passband_gain, not both
If norm_factor is not specified, it is automatically calculated so that
    A_0 * abs(Hp(s)) = 1 at f0
as per the SEED recommendation

If A_0 is automatically calculated, it will change if the user changes ref_freq
The user can also change the input_units after creating the class, in which
case the poles, zeros, AUTOMATIC norm_factor and gain will all
be changed if the starting and ending input_units are in self.known_units().
The norm_factor will NOT be changed if it was manually entered.
"""
import numpy as np
import json
import copy

from obspy.core.inventory import read_inventory
from obspy.core.inventory.response import PolesZerosResponseStage
from .logger import log


class PAZ():
    """
    Object-based implementation of Pole-Zero representation

    For representing a full instrument for simulations, etc
    """
    converter = {"m":      [1.,    2],
                 "m/s":    [1.,    1],
                 "m/s^2":  [1.,    0],
                 "mm":     [1.e-3, 2],
                 "mm/s":   [1.e-3, 1],
                 "mm/s^2": [1.e-3, 0],
                 "nm":     [1.e-9, 2],
                 "nm/s":   [1.e-9, 1],
                 "nm/s^2": [1.e-9, 0],
                 "gal":    [1.e-2, 0],   # 1 cm/s^2
                 "mgal":   [1.e-5, 0],
                 "ugal":   [1.e-8, 0]}

    def __init__(self, poles=[], zeros=[], passband_gain=None,
                 input_units='m/s', output_units='counts', ref_freq=1.,
                 gain=None, norm_factor=None):
        """
        Can either specify gain or passband_gain, not both

        Frequency units are radians

        :param poles: poles (rad/s)
        :param zeros: zeros (rad/s)
        :param passband_gain: gain in the instrument passband
        :param ref_freq: reference frequency (Hz) for passband gain
        :param input_units: units (usually physical) at the sensor entry
        :param output_units: units as stored
        :param gain: pole-zero gain (NOT passband gain)
        :param norm_fact: A0 normalization factor (blocks automatic
            calculation)
        """
        assert gain is None or passband_gain is None,\
            "both gain and passband_gain cannot be declared"
        if gain is None and passband_gain is None:
            passband_gain = 1.
        self.poles = np.array(poles)
        self.zeros = np.array(zeros)
        self._input_units = input_units
        self.output_units = output_units
        self.ref_freq = ref_freq
        self._norm_factor = norm_factor
        if gain is not None:
            self.gain = gain
        else:
            self.gain = passband_gain * self.norm_factor

    @property
    def passband_gain(self):
        return self.gain / self.norm_factor

    @property
    def norm_factor(self):
        """ A0 normalization factor """
        if self._norm_factor is not None:
            return self._norm_factor
        else:
            return self.calc_norm_factor()

    @norm_factor.setter
    def norm_factor(self, value):
        self._norm_factor = value

    @property
    def input_units(self):
        return self._input_units

    @input_units.setter
    def input_units(self, value):
        """
        Convert input (physical) units to a different value

        If newunits and self.input_units are one of the
        DEFINED UNITS (see _ref_meter()), changes gain, poles and zeros.
        """
        assert isinstance(value, str), 'value is not a character string'
        self._unit_conversion(value)
        self._input_units = value

    def calc_norm_factor(self):
        """
        Calculate A0 normalization factor

        Returns A0 such that A0 * prod(s-z_n)/prod(s-p_n) = 1 at ref_freq
        (SEED recommendation)
        """
        s = 2 * np.pi * self.ref_freq * 1j
        val = 1
        for p in self.poles:
            val *= (s - p)
        for z in self.zeros:
            val /= (s - z)
        val = abs(val)
        return val

    def __eq__(self, obj, places=2):
        """
        Return true if objects are equal to within places decimal places

        :param places: decimal place precision used for gain
        """
        diff = abs(self.gain/obj.gain - 1) * 10**places
        if diff > 1:
            return False
        if not len(self.zeros) == len(obj.zeros):
            return False
        if not len(self.poles) == len(obj.poles):
            return False
        for z1, z2 in zip(self.zeros, obj.zeros):
            if not z1 == z2:
                return False
        for p1, p2 in zip(self.poles, obj.poles):
            if not p1 == p2:
                return False
        if not self.input_units.upper() == obj.input_units.upper():
            return False
        if not self.output_units.upper() == obj.output_units.upper():
            return False
        if not self.ref_freq == obj.ref_freq:
            return False
        return True

    def __str__(self):
        s = 'PAZ: '
        s += f'gain={self.gain:.4g}, '
        s += f'poles={self.poles}, '
        s += f'zeros={self.zeros}, '
        s += f'[passband_gain={self.passband_gain:.3g} '
        s += f'{self.output_units}/{self.input_units} '
        s += f'at {self.ref_freq:.3g}Hz, '
        s += f'A0={self.norm_factor:.4g}]'
        return s

    def copy(self):
        return copy.deepcopy(self)

    def to_obspy(self):
        """
        Return obspy paz dict

        obspy is not clear about what is in the "gain" and "sensitivity"
        keys.  I tried self.norm_factor and self.gain, respectively, but
        this makes their product depend on f_ref.  Switched to:
            self.norm_factor, self.passband_gain
        """
        return {'poles': list(self.poles), 'zeros': list(self.zeros),
                'f_ref': self.ref_freq, 'gain': self.norm_factor,
                'sensitivity': self.passband_gain}

    def _unit_conversion(self, new_input_units):
        """
        Convert from self.input_units to new_input_units
        """
        igain, inz = self._ref_meter(self.input_units)
        ogain, onz = self._ref_meter(new_input_units)

        if igain is None or ogain is None:
            log('Did not convert', 'warning')
            return

        self.gain *= ogain / igain
        relnz = onz - inz
        if relnz > 0:
            self.zeros = np.append(self.zeros, np.zeros(relnz))
        elif relnz < 0:
            relnp = -relnz
            i_zero_zeros = np.where(np.absolute(self.zeros) == 0)[0]
            if len(i_zero_zeros) > 0:
                if len(i_zero_zeros) >= relnp:
                    self.zeros = np.delete(self.zeros, i_zero_zeros[:relnp])
                    relnp = 0
                else:
                    self.zeros = np.delete(self.zeros, i_zero_zeros)
                    relnp -= len(i_zero_zeros)
            self.poles = np.append(self.poles, np.zeros(relnp))

    @classmethod
    def read_JSON_PZ(cls, filename):
        """
        read JSON paz format
        """
        with open(filename) as f:
            paz = json.load(f)
        poles, zeros = [], []
        for p in paz['poles']:
            poles.append(float(p[0]) + float(p[1])*1.j)
        for z in paz['zeros']:
            zeros.append(float(z[0]) + float(z[1])*1.j)
        val = cls(poles=poles, zeros=zeros,
                  passband_gain=1./float(paz['sensitivity']),
                  ref_freq=float(paz['f_ref']),
                  input_units='nm',
                  output_units='counts')
        return val

    @classmethod
    def read_STATIONXML(cls, filename, component):
        """
        read JSON paz format
        """
        inv = read_inventory(filename, 'STATIONXML')
        resp = inv.select(channel='*' + component)[0][0][0].response
        s_pz = []
        s_others = []
        stages = resp.response_stages
        sens = resp.instrument_sensitivity
        for s in stages:
            if isinstance(s, PolesZerosResponseStage):
                s_pz.append(s)
            else:
                s_others.append(s)
        if len(s_pz) > 1:
            paz_base = s_pz[0]
            for s in s_pz[1:]:
                paz_base.poles.extend(s.poles)
                paz_base.zeros.extend(s.zeros)
                paz_base.normalization_factor *= s.normalization_factor
            resp.response_stages = [paz_base] + s_others
        paz = resp.get_paz()
        val = cls(poles=paz.poles, zeros=paz.zeros,
                  passband_gain=sens.value,
                  ref_freq=sens.frequency,
                  input_units=sens.input_units,
                  output_units=sens.output_units)
        return val

    @classmethod
    def read_SACPZ(cls, filename):
        """
        read SAC PZ format

        input units are assumed to be meters (does not look at comments)
        The file-reading code is copied from obspy.io.sac.sacpz.attach_paz()
        """
        with open(filename) as f:
            poles = []
            zeros = []
            while True:
                line = f.readline()
                if not line:
                    break
                # lines starting with * are comments
                if line.startswith('*'):
                    continue
                if line.find('ZEROS') != -1:
                    a = line.split()
                    noz = int(a[1])
                    for _k in range(noz):
                        line = f.readline()
                        a = line.split()
                        if line.find('POLES') != -1 or \
                           line.find('CONSTANT') != -1 or \
                           line.startswith('*') or not line:
                            while len(zeros) < noz:
                                zeros.append(complex(0, 0j))
                            break
                        else:
                            zeros.append(complex(float(a[0]), float(a[1])))

                if line.find('POLES') != -1:
                    a = line.split()
                    nop = int(a[1])
                    for _k in range(nop):
                        line = f.readline()
                        a = line.split()
                        if line.find('CONSTANT') != -1 or \
                           line.find('ZEROS') != -1 or \
                           line.startswith('*') or not line:
                            while len(poles) < nop:
                                poles.append(complex(0, 0j))
                            break
                        else:
                            poles.append(complex(float(a[0]), float(a[1])))
                if line.find('CONSTANT') != -1:
                    a = line.split()
                    constant = float(a[1])
        val = cls(poles=poles, zeros=zeros,
                  gain=constant,
                  input_units='m',
                  output_units='counts')
        return val

    @classmethod
    def read_GSE(cls, filename):
        """
        read a GSE response file

       From specification in GSE_2.0.pdf (Group of Scientific Expernts, 1995)
        :param name:    Name of GSE2 response file
        :returns: Structure with all parameters
        """
        # read_GSE.m:6
        info = dict()
        with open(filename, 'r') as fic:
            for line in fic:
                # CAL2 line is specified in Table 13, page 69
                if line[:4] == 'CAL2':
                    info['station'] = line[5:10]
                    info['type'] = line[20:26].strip()
                    info['sensitivity'] = float(line[27:37])
                    info['f_ref'] = 1. / float(line[38:45])
                    info['sampling_rate'] = float(line[46:56])
                # PAZ2 line is specified in Table 14, page 70
                elif line[:4] == 'PAZ2':
                    info['norm_const'] = float(line[10:25])
                    num_poles = float(line[40:43])
                    num_zeros = float(line[44:47])
                    info['poles'] = []
                    info['zeros'] = []
                    for p in np.arange(num_poles):
                        A = next(fic).split()
                        info['poles'].append(np.complex(float(A[0]),
                                                        float(A[1])))
                    for z in np.arange(num_zeros):
                        A = next(fic).split()
                        info['zeros'].append(np.complex(float(A[0]),
                                                        float(A[1])))
                # DIG2 line is specified in Table 17, page 72
                elif line[:4] == 'DIG2':
                    info['amplifier_gain'] = float(line[8:23])
                    info['digitizer'] = line[36:].strip()
                else:
                    continue
        info['f_ref'] = info.get('f_ref', 1.)
        val = cls(poles=info['poles'], zeros=info['zeros'],
                  passband_gain=1./info['sensitivity'],
                  ref_freq=info['f_ref'],
                  input_units='nm',
                  output_units='counts')
        # if info['A0'] / val.norm_factor < 0.99\
        #         or info['A0']/val.norm_factor > 1.01:
        #     log('self.norm_factor is different from given: {:.3g} vs {:.3g}'
        #         .format(val.norm_factor, info['A0']), 'warning')
        return val

    @classmethod
    def read_Baillard_PZ(cls, filename):
        """
        read Baillard PZ format

        A numbers-only version of the GSE format
        """
        info = {}
        with open(filename, 'r') as foc:
            i = 0
            a, info['poles'], info['zeros'] = [], [], []
            for tline in foc:
                if i < 4:
                    a.append(float(tline))
                elif i == 4:
                    n_poles = float(tline)
                elif i == 5:
                    n_zeros = float(tline)
                elif i >= 6 and i < (6 + n_poles):
                    A = tline.split()
                    info['poles'].append(np.complex(float(A[0]), float(A[1])))
                elif i >= 6 + n_poles and i < 6 + n_poles + n_zeros:
                    A = tline.split()
                    info['zeros'].append(np.complex(float(A[0]), float(A[1])))
                else:
                    info['amplifier_gain'] = float(tline)
                i += 1
        info['sensitivity'] = a[0]
        info['f_ref'] = 1 / a[1]
        info['sampling_rate'] = a[2]
        info['A0'] = a[3]
        val = cls(poles=info['poles'], zeros=info['zeros'],
                  passband_gain=1/info['sensitivity'],
                  ref_freq=info['f_ref'],
                  # input_units='m/s',
                  input_units='nm',
                  output_units='counts')
        # if info['A0'] / val.norm_factor < 0.99\
        #         or info['A0'] / val.norm_factor > 1.01:
        #     log('self.norm_factor is different from given: {:.3g} vs {:.3g}'
        #         .format(val.norm_factor, info['A0']), 'warning')
        return val

    def known_units(self):
        """
        Return known input_units (for automatic conversions
        """
        return list(self.converter.keys())

    def _ref_meter(self, units):
        """
        Return conversions needed to convert to meters

        :returns: gain_factor, number of zeros
        """
        if units.lower() in self.converter:
            return (self.converter[units.lower()][0],
                    self.converter[units.lower()][1])
        else:
            log('Unknown unit for conversion: {}, known values are "{}"'
                .format(units, ', '.join(self.converter.keys())))
            return None, None

# function ps=butterPoles(nPoles,om)
#     if length(nPoles)>1, error('length(nPoles) > 1'); end
#     if nPoles==0,ps=[]; return; end
#     if length(om)>1, error('length(om) > 1'); end
#     phase=((2*(1:nPoles)+nPoles-1)*pi)/(2*nPoles);
#     normPoles=exp(1i*phase);
#     %phase90=phase*180/pi   % [phase in degrees for debugging
#     ps=om*normPoles;
#     ps=ps(:)';  % Force row matrix
# end
#
# classdef ts_response
#     % Class for instrument responses
#     %
#     % RESPONSE
#     %   The response specifies how your data was transformed from
#     %   physical ("input") values to the stored ("output") values:
#     %
#     %   outUnits = inUnits * resp, where:
#     %
#     %                    (s-zero1)*(s-zero2)*(s-zero3)*...
#     %	resp(s) = pzgain*----------------------------------
#     %                    (s-pole1)*(s-pole2)*(s-pole3)*...
#     %
#     % 	where s=1i*2*pi*frequency
#     %
#     %   The reponse is the transfer function FROM the input units TO
#     %   the output units.
#     %   We calculate spectra and convert time series to input units by
#     %   dividing the stored data by this response.
#     %
#     % GAIN VERSUS "PASSBAND" GAIN
#     %   The pole-zero gain stored in this class is not the "passband gain"
#     %   given for standard lowpass, highpass or bandpass filters, because
#     %   transfer functions do not necessarily have a flat passband in which
#     %   this gain makes sense.
#     %   If you want to specify the passband gain, create the object
#     %     1) using pole-zero notation and specifying omega_ref, or
#     %     2) using Butterworth filter paramters
#     %   for details, see the help for the creator method.
#     %
#     % ts_response properties:
#     %   gain     - Instrument pole-zero gain
#     %   poles    - Instrument poles (radians)
#     %   zeros    - Instrument zeros (radians)
#     %   inUnits  - Units of data at the instrument's input
#     %              (e.g., physical units)
#     %   outUnits - Units of stored data (e.g., digital units)
#     %
#     % ts_response methods:
#     %   Creator
#     %      ts_response    - Creator
#     %   Miscellaneous
#     %      getNorm        - Return normalization factor and gain at a freq
#     %      plot           - Plot the response
#     %      response       - Return response at given frequencies
#     %   Reading from standard format
#     %      Use the ts_response creator with the filetype as the first
#     %      argument
#     %   Writing to standard formats
#     %      writeGSE2_PZ   - Write the response in GSE pole-zero format (v2.0)
#     %      writeSAC_PZ    - Write the response to a SAC pole-zero file
#     %      writeSEISAN_PZ - Write the response to a SEISAN pole-zero file
#     %   Overloaded Matlab methods
#     %      eq             - Implements "="
#     %      ne             - Implements "~="
#     %      mtimes         - Implements "*"
#     %
#     % EXAMPLES:
#     %   % Bandpass filter with 2 poles at 4.5 Hz and 2 at 100 Hz
#     %   tsr=ts_reponse('BW','m/s^2','du',1000,[100 100],[4.5 4.5]);
#     %   plot(tsr);
#     %
#     %   % The same thing in pole-zero notation
#     %   tsr=ts_response('PZ','m/s^2','du',3.94784e8,-2*pi*[4.5 4.5 100 100],
#     %                   [0 0]);
#     %   plot(tsr)
#     %
#     %   % Read in from a SAC pole-zero file
#     %   tsr=ts_response('SACPZ','','','MYCHANNEL.SACPZ');
#     %   plot(tsr)
#     %
#     %   % ALMOST the same thing in pole-zero notation, specifying gain at
#     %   %  a frequency in the passband (but the passband isn't perfectly
#     %   %  flat, so the response isn't exactly the same
#     %   tsr=ts_response('PZ','m/s^2','du',1000,-2*pi*[4.5 4.5 100 100],
#     %                   [0 0],2*pi*20);
#     %   plot(tsr)
#
#     %% PROPERTIES
#     properties
#         gain=1;              % Instrument pole-zero gain (NOT passband gain)
#         poles=[];            % Poles
#         inUnits='';         % Units of data that you want to see
#         outUnits='';        % Units of data as they are stored
#     end
#     properties (Dependent=true)
#         zeros               % zeros
#     end
#
#     properties (Access=private, Hidden=true)
#         privateZeros=[];           % True storage of zeros (necessary so that
#         % matlab doesn't confuse with the function zeros())
#     end
#
#     %% PUBLIC METHODS
#     methods
#         %% CREATOR
#         function [obj,extra]=ts_response(inpType,iu,ou,argA,argB,argC,
#                                          argD,argE)
#             % Create a ts_response object
#             %
#             % Usage:
#             %   obj= ts_response();
#             %   obj= ts_response('SACPZ',inUnits,outUnits,filename);
#             %   obj= ts_response('GSE2PZ',filename);
#             %   obj= ts_response('PZ',inUnits,outUnits,gain,poles,zeros);
#             %   obj= ts_response('PZ',inUnits,outUnits,gain,poles,zeros,
#             %                    om_ref); % ""    ""   ""
#             %   obj= ts_response('BW',inUnits,outUnits,gain,nHP,hpFrq,nLP
#             %                    ,lpFrq);% Butterworth Filter
#             %
#             % Inputs:
#             %   type	Input type:
#             %               'SACPZ': SAC Pole-Zero file        %
#             %               'GSE2PZ': GSE2 Pole-Zero file
#             %               'PZ':    Pole-Zero Notation
#             %               'BW':    Butterworth Filter
#             %               'Units': Unit Conversion
#             %   inUnits	The input units of the data, that is,
#             %                   the units that you want spectra, etc to be
#             %                   calculated in
#             %   outUnits The units in which the data are stored;
#             %
#             %  SAC POLE-ZERO FILE
#             %   filename        Name of a SAC pole-zero file (poles & zeros
#             %                   in rad/s)
#             %   inUnits         Automatically set to 'nm' if empty
#             %  GSE2 POLE-ZERO FILE
#             %   filename        Name of a GSE2 pole-zero file (poles & zeros
#             %                   in rad/s)
#             %  POLE-ZERO NOTATION
#             %   gain            Pole-zero gain (NOT the passband gain,
#             %                   unless omega_ref is specified)
#             %   poles           Pole-zero poles (radians/s)
#             %   zeros           Pole-zero zeros (radians/s)
#             %   om_ref          Frequency (radians/s) at which  gain applies
#             %                    (useful for implementing a "passband gain")
#             %  BUTTERWORTH FILTER
#             %   gain            Passband gain
#             %   nLP             Number of highpass poles
#             %   lpFrqs          Highpass corner frequency (Hz)
#             %   nLP             Number of lowpass poles
#             %   lpFrqs          Lowpass corner frequency (Hz)
#             %
#             % Output:
#             %   obj           the ts_response object
#             %
#             % Following convention, poles and zeros are expressed in radians
#             %         and Butterworth filter frequencies are in Hz.  Here are
#             %         some conversion factors.
#             %            radians = 2*pi*Hz = 2*pi/Period = 1/RC
#             %            Hz = radians/2pi = 1/Period = 1/(2*pi(RC)
#             %         where Period is the instrument "natural period" and RC
#             %         is the resistance*capicitance for each stage of an
#             %         analog electric filter.
#             %
#             %
#             % Examples:
#             %	resp = ts_response('Units','m/s^2','nm');
#             %	resp = ts_response('SACPZ','','','myPZs');
#             %	resp = ts_response('BW','Pa','V',200,1,50);
#             %% Check input arguments
#             ts_debug(1,'Starting with %d arguments\n',nargin);
#             extra=[];
#             if nargin==0;  return;    % 0 arguments => Create empty object
#             % 1 argument: ts_response object => copy
#             elseif ((nargin==1) && isa(inpType,'ts_response'))
#                 obj=inpType; return;
#             elseif ~ischar(inpType),
#                 error(' The first argument is not a string (type ''doc '
#                       'ts_response'')');
#             end
#             if nargin >=4,
#                 obj.inUnits = iu;
#                 obj.outUnits = ou;
#             end
#
#             %% Create poles, zeros and gain according to the input format
#             normalizer=1;
#             switch upper(inpType),
#                 case 'SACPZ'
#                     if nargin ~= 4
#                         error('ts_response(''SACPZ'',...) needs 4 arguments,'
#                               ' you gave %d',nargin)
#                     end
#                     [g,ps,zs] = readSACPZfile(argA);
#                     if isempty(obj.inUnits)
#                         obj.inUnits='nm';    % Standard for SAC PZ files
#                     end
#                 case 'GSE2PZ'
#                     if nargin ~= 1 && nargin ~=2,
#                         error('ts_response(''GSE2PZ'',...) needs 1 or 2 '
#                               'arguments, you gave %d',nargin)
#                     end
#                     if ~exist('iu','var'); argA=[];
#                     else argA=iu;
#                     end
#                     [g,ps,zs,~] = readGSE2PZfile(argA);
#
#                     obj.inUnits='nm'; obj.outUnits='counts';
#                 case 'PZ',
#                     if nargin~=6 && nargin~=7
#                         error('ts_response(''PZ'',...) needs 6 or 7 '
#                               'arguments, you gave %d',nargin)
#                     end
#                     g=argA;
#                     ps=argB(:)';
#                     zs=argC(:)';
#                     if exist('argD','var')
#                         omega_norm = argD;
#                         s=1i*omega_norm;
#                         for p = ps, normalizer = normalizer*(s-p); end
#                         for z = zs, normalizer = normalizer/(s-z); end
#                         normalizer=abs(normalizer);
#                         ts_verbose(1,'Calculated normalizer = %.2f at '
#                                    '%g radians',abs(normalizer),omega_norm);
#                         g=g*normalizer;
#                     end
#                 case {'BADBW'}, % Not a real butterworth, all poles real
#                     if nargin~=6
#                         error('ts_response(''badBW'',...) needs 6 arguments,'
#                               ' you gave %d',nargin)
#                     end
#                     hpOmega=2*pi*argC;
#                     lpOmega=2*pi*argB;
#                     nLP = length(lpOmega);
#                     nHP = length(hpOmega);
#                     ps=zeros(1,nLP+nHP); zs=zeros(1,nHP);
#                     for i = 1:length(hpOmega),
#                         ps(i) = -hpOmega(i);
#                         zs(i) = 0;
#                     end
#                     for i=1:length(lpOmega),
#                         ps(i+nHP) = -lpOmega(i);
#                         normalizer = normalizer*lpOmega(i);
#                     end
#                     normalizer=abs(normalizer);
#                     g = argA*normalizer;
#                 case {'BW','BUTTERWORTH'},
#                     assert(nargin==8,
#                            'ts_response(''BW'',...) needs 8 arguments, you '
#                             'gave %d',nargin);
#                     nHP=argB;
#                     assert(length(nHP)==1,'length(nHP)~=1');
#                     assert(isnumeric(nHP),'nHP is not numeric');
#                     hpOmega=2*pi*argC;
#                     if nHP==0,
#                             assert(isempty(hpOmega),
#                                    '~isempty(hpOmega) for nHP==0');
#                     else assert(length(hpOmega)==1,'length(hpOmega)~=1'); end
#                     assert(isnumeric(hpOmega),'hpOmega is not numeric');
#                     nLP=argD;
#                     assert(length(nLP)==1,'length(nLP)~=1');
#                     assert(isnumeric(nLP),'nLP is not numeric');
#                     lpOmega=2*pi*argE;
#                     if nLP==0,
#                         assert(isempty(lpOmega),
#                                '~isempty(lpOmega) for nLP==0');
#                     else assert(length(lpOmega)==1,'length(lpOmega)~=1');end
#                     assert(isnumeric(lpOmega),'lpOmega is not numeric');
#
#                     ps=[butterPoles(nHP,hpOmega) butterPoles(nLP,lpOmega)];
#                     zs=zeros(1,nHP);
#                     if ~isempty(lpOmega)
#                         normalizer=abs(lpOmega^nLP);
#                     else normalizer=1;
#                     end
#                     g = argA*normalizer;
#                 otherwise,
#                     error('Invalid type: %s.  Valid types are ''PZ'', '
#                           '''SACPZ'', ''BW'' and ''UNITS''',inpType);
#             end
#
#             %% Save poles, zeros and gain
#             obj.gain = g;
#             if length(obj.gain) ~= 1,
#                 error('gain must be a scalar');
#             end
#             obj.poles = ps(:)';
#             obj.privateZeros = zs(:)';
#             ts_debug(1,'created ts_response object:\n%s',char(obj));
#
#         end    % Creator function ts_response
#
#         %% SET METHODS (ELIMINATE BAD INPUTS)
#         function obj = set.gain(obj,value)
#             if ~isa(value,'numeric'), error('value must be numeric'); end
#             obj.gain=value;
#         end
#
#         function obj = set.poles(obj,value)
#             if ~isa(value,'numeric'), error('value must be numeric'); end
#             [M,N]=size(value);
#             if M~=1 && N~=1, error('value must be one-dimensional'); end
#             obj.poles=value;
#         end
#
#         function obj = set.zeros(obj,value)
#             if ~isa(value,'numeric'), error('value must be numeric'); end
#             [M,N]=size(value);
#             if M~=1 && N~=1, error('value must be one-dimensional'); end
#             obj.privateZeros=value;
#         end
#
#         function obj = set.outUnits(obj,value)
#             if ~isa(value,'char'),
#                 error('value must be a character string'); end
#             obj.outUnits=value;
#         end
#
#         %% OTHER METHODS
#         function result = eq(a,b)
#             % Implements a==b
#             % returns 1 if a is the same as b, 0 otherwise
#             if (~isa(b,'ts_response')),
#                 error('second object must also be ts_response');
#             end
#             result=1;
#
#             if a.gain ~= b.gain,					result=0; return; end
#             if length(a.poles) ~= length(b.poles),	result=0; return; end
#             if length(a.zeros) ~= length(b.zeros),	result=0; return; end
#             for i=1:length(a.poles),
#                 if a.poles(i) ~= b.poles(i),		result=0; return; end
#             end
#             for i=1:length(a.zeros),
#                 if a.zeros(i) ~= b.zeros(i),		result=0; return; end
#             end
#         end
#         function result = ne(a,b)
#             % Implements obja ~= objb
#             % returns 1 if a not the same as b, 0 otherwise
#
#             result = ~(a==b);
#         end
#         plot(varargin)                    % DEFINED IN SEPARATE FILE
#
#         function fpout = mtimes(fpa,fpb)
#             % Implements obja * objb
#             %
#             % Example: newfilt=filta*filtb;
#             %
#             % Multiplies the gains, concatenates the poles and zeros
#             %
#             if isa(fpa,'ts_response') && isa(fpb,'ts_response'),
#                 fpout=ts_response();
#                 fpout.gain=fpa.gain*fpb.gain;
#                 fpout.poles=[fpa.poles fpb.poles];
#                 fpout.privateZeros=[fpa.privateZeros fpb.privateZeros];
#                 if strcmp(fpa.inUnits,fpb.outUnits)
#                     fpout.inUnits=fpb.inUnits;
#                     fpout.outUnits=fpa.outUnits;
#                 elseif strcmp(fpb.inUnits,fpa.outUnits)
#                     fpout.inUnits=fpa.inUnits;
#                     fpout.outUnits=fpb.outUnits;
#                 % One filter does not not change units, use units from other
#                 elseif isempty(fpa.inUnits)&&isempty(fpa.outUnits)
#                     fpout.inUnits=fpb.inUnits;
#                     fpout.outUnits=fpb.outUnits;
#                 % One filter does not not change units, use units from other
#                 elseif isempty(fpb.inUnits)&&isempty(fpb.outUnits)
#                     fpout.inUnits=fpa.inUnits;
#                     fpout.outUnits=fpa.outUnits;
#                 else
#                     ts_warning('One response''s input Units does not match '
#                                'the other response''s output Units.  Setting'
#                                ' input and output Units to UNKNOWN');
#                     fpout.inUnits='UNKNOWN';
#                     fpout.outUnits='UNKNOWN';
#                 end
#             elseif isa(fpa,'double') && length(fpa)==1,
#                 fpout=fpb; mulfact=fpa;
#                 fpout.gain=fpout.gain*mulfact;
#             elseif isa(fpb,'double') && length(fpb)==1,
#                 fpout=fpa; mulfact=fpb;
#                 fpout.gain=fpout.gain*mulfact;
#             else
#                 error('Both objects must be ts_response objects, or one a '
#                       'scalar');
#             end
#         end
#         function resp = response(obj,fvect)
#             % Return instrument response for a range of frequencies
#             %
#             % Usage:
#             %	resp = response(filtp,freqs);
#             %	Calculates response of a set of filter parameters
#             %		filtparm is an ts_response object
#             %		freqs is a vector of frequencies (1/s)
#             if ~exist('fvect','var'),
#                 maxf=100;fstep=maxf/100; fvect=fstep:fstep:maxf;
#             end
#             onevect = ones(size(fvect)); omega = 2*pi*fvect;
#             i = sqrt(-1);
#             % Gain
#             resp = obj.gain * onevect;
#             % Poles
#             for p = obj.poles,
#                 resp = resp ./ (i*omega - p);
#             end
#             % Zeros
#             for z = obj.privateZeros,
#                 resp = resp .* (i*omega - z);
#             end
#         end
#     end % public methods
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function [gain,ps,zs] = readSACPZfile(fname)
#     % read in pole-zero information from a SAC pole-zero file:
#     % Frequencies are specified in omega format
#     % Poles are specified by the line POLES n, where n is the # of poles,
#     % followed by m lines with the real and imaginary parts of each pole.
#     % If m<n, remaining poles are assumed to be at the origin
#     % Zeros are specified in the same way using the line ZEROS n
#     % A multiplicative constant to add is specified using CONSTANT x
#     % Up to 15 poles and 15 zeros, I THINK SAC convention is that the resp
#     % is in velocity (m/s).  AnyPlot doesn't assume anything, you  enter the
#     % units as text in the "chCorrectedUnits" variable
#
#     t.gain = 1;
#     t.poles = [];
#     t.privateZeros = [];
#
#     inPZ=0;	% set to 1 if reading poles, set to 2 if reading zeros
#     nPZ=0; iPZ=0;
#     fid = fopen(fname,'r');
#     if fid==0, error('Pole-zero file %s not found',fname); end
#     ts_debug(1,'Reading pole-zero file %s\n',fname);
#     while 1,
#         tline=fgetl(fid);
#         if ~ischar(tline), break, end
#         if length(deblank(tline))<=1, break,end
#         ts_debug(1,tline);
#         if tline(1)=='*',        % comment line
#             continue
#         end
#         A = strread(tline,'%s');
#         %A{1}
#         %A = textscan(tline,'%s')
#         %A{1}
#         %ts_debug(1,'A=%s\n',celldisp(A));
#         switch A{1},
#             case 'POLES',
#                 t = cleanUpOld(t,inPZ,nPZ,iPZ);
#                 inPZ = 1;
#                 nPZ = str2double(A{2});
#                 iPZ = 0;
#             case 'ZEROS',
#                 t = cleanUpOld(t,inPZ,nPZ,iPZ);
#                 inPZ = 2;
#                 nPZ = str2double(A{2});
#                 iPZ = 0;
#             case 'CONSTANT'
#                 t = cleanUpOld(t,inPZ,nPZ,iPZ);
#                 t.gain = sscanf(tline,'CONSTANT %f');
#             case '*'
#                 t = cleanUpOld(t,inPZ,nPZ,iPZ);
#                 t.gain = sscanf(tline,'CONSTANT %f');
#             otherwise;
#                 omega = str2double(A{1}) + sqrt(-1)*str2double(A{2});
#                 % f = omega/(2*pi);
#                 if inPZ==1,
#                     t.poles=[t.poles omega];
#                     iPZ=iPZ+1;
#                 elseif inPZ==2,
#                     t.privateZeros=[t.privateZeros omega];
#                     iPZ=iPZ+1;
#                 else
#                     error('non-command line occured while not in list mode')
#                 end
#         end
#     end
#
#     ps = t.poles(:)';
#     zs = t.privateZeros(:)';
#     gain=t.gain;
# end
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function [gain,ps,zs,GSEprops] = readGSE2PZfile(fname)
#     % read response from a GSE 2.0 pole-zero file
#     % Usage: [gain,ps,zs,GSEprops] = readGSE2PZfile(fname)
#     %
#     %   gain,ps,zs are the gain, poles and zeros
#     %   GSEprops is a structure with the following fields:
#     %           startDate   Date from which these values apply
#     %           refFreq     Freq at which to give the system response (nm/c)
#     %           station     Station name
#     %           comp        Component name (up to three characters)
#     %           inSPS       Sample rate at digitizer input
#     %           outSPS      Sample rate at digitizer output
#     %           sensorType  Sensor type (i.e., 'STS2' 'T240' 'L-28')
#     %           digType     Digitizer type (i.e. 'DM-24', 'L-CHEAPO')
#     %           digGain     Digitizer gain (counts/V, object's gain will be
#     %                       divided by this
#     %                           to calculate the sensor gain)
#     % The GSE 2.0 response format is:
#     % "CAL2" station comp sensorType refgain(nm/c) refFreq(Hz) outputSPS
#     %        startDate
#     % "PAZ2  1 V"  sensorNormConst(V/nm) nPoles  nZeros  "Laplace transform"
#     % real(Pole1)  imag(Pole1)
#     % real(Pole2)  imag(Pole2)
#     % ....
#     % real(Zero1)  imag(Zero1)
#     % read(Zero2)  imag(Zero2)
#     % ...
#     % "DIG2  2" digitizerGain(c/V)  inputSPS  digitizerType
#     %
#     % Frequencies are specified in radians
#     %
#     % The file specifications come from
#     % http://www.seismo.ethz.ch/prod/autodrm/manual/CRP243.OpsAnx3.A4.pdf
#     %
#     % Note that the GSE2.1 format is different enough to be unreadble by a
#     % 2.0 reader, but has no tag to distinguish from 2.0!
#     % http://www.seismo.ethz.ch/prod/autodrm/manual/provisional_GSE2.1.pdf
#
#     if ~exist('fname','var')
#         fname=[]; end
#     if isempty(fname);
#         [FileName,PathName]=uigetfile({'*_GSE','GSE PZ files'},
#                                       'Select a GSE 2.0 Pole-Zero File');
#         fname=fullfile(PathName,FileName);
#     end
#     fid=fopen(fname,'r');
#     if fid==0, error('GSE2.0 PZ file %s not found',fname);
#     else
#         ts_debug(1,'Reading pole-zero file %s\n',fname);
#     end
#
#     GSEprops.comments='';
#     line=fgetl(fid);
#     % LINE 1: GENERAL INFORMATION (none of this is stored, but total gain
#     % is compared with component gains
#     if ~strcmp(line(1:4),'CAL2'),
#         error('Line 1 does not start with CAL2');
#     else
#         station=strtrim(line(6:10));
#         comp=   strtrim(line(12:14));
#         %auxID=  strtrim(line(16:19));
#         instType=strtrim(line(21:26));
#         calib=   str2double(line(28:37));  % system sens (nm/count) @ ref pd
#         %calib=   str2double(line(28:42));    % GSE 2.1
#         calPer=  str2double(line(39:45));    % calibration reference pd (s)
#         %calPer=  str2double(line(44:50));    % GSE2.1
#         outSPS=  str2double(line(47:56));
#         %outSPS=  str2double(line(52:62));    % GSE 2.1
#         onDate=  strtrim(line(58:67));
#         %onDate=  strtrim(line(64:73));    % GSE 2.1
#         onTime=  strtrim(line(69:73));
#         %onTime=  strtrim(line(75:79));    % GSE 2.1
#         %offDate= strtrim(line(75:84));
#         %offDate= strtrim(line(81:90));    % GSE 2.1
#         %offTime= strtrim(line(86:90));
#         %offTime= strtrim(line(92:96));    % GSE 2.1
#
#         ts_verbose(0,'CAL2 line: station comp instType refGain(nm/cnt) @Hz '
#                    '  outSPS   startDate');
#         ts_verbose(0,'  %-5s   %-3s  %-7s  %-12g    %-5g %-8g %s %s', ...
#             station,comp,instType,calib,1/calPer,outSPS,onDate,onTime);
#     end
#
#     %% READ IN RESPONSE INFORMATION: earth input is assumed to be in nm
#     %% displacement
#     line=fgetl(fid);
#     pazDone=0; digDone=0;
#     while 1,
#         if isempty(line), continue
#         else
#             c=strtrim(line); % remove white spaces
#             if c(1)=='('     % Comment line
#                 if c(end)~=')'
#                     error('Comment not closed by parentheses: "%s"',line);
#                 end
#                 ts_verbose(0,'COMMENT: %s',c);
#                 comment=strtrim(c(2:end-1));
#                 if isempty(GSEprops.comment)
#                     GSEprops.comment=comment;
#                 else
#                     if ~iscell(GSEprops.comment)
#                         GSEprops.comment={GSEprops.comment};
#                     end
#                     GSEprops{length(GSEprops.comment)+1}=comment;
#                 end
#                 continue
#             end
#         end
#         switch(line(1:4))
#             case 'FIR2', error('Finite Impulse response not implemented')
#             case 'GEN2', error('Generic Response not implemented')
#             case 'FAP2',
#                 error('Frequency, Amplitude, Phase Response not implemented')
#             case 'PAZ2'     % POLE-ZERO RESPONSE (usually SEISMOMETER)
#                             % [Same for GSE 2.0 and 2.1]
#                 if pazDone,
#                     error('Code can't handle Multiple Pole-Zero responses');
#                 end
#                 %sNum=str2double(line(6:7));    % Stage sequence number
#                 oUnits=line(9);
#                 if oUnits ~= 'V',
#                    error('outUnits are %s, function only written for "V"'
#                          ,oUnits); end
#                 sFactor=str2double(line(11:25));
#                 %deci=str2double(line(27:30));      % Unused
#                 %gCorr=str2double(line(32:39));   % Group correction applied
#                                                   % seconds)
#                 nPoles=str2double(line(41:43));
#                 nZeros=str2double(line(45:47));
#                 descript=strtrim(line(49:end));
#                 % sFactor = ounits/nm
#                 j=sqrt(-1);
#                 ps=zeros(1,nPoles);
#                 for i = 1:nPoles,
#                     line=fgetl(fid);
#                     val=sscanf(line,'%f %f');
#                     ps(i)=val(1)+j*val(2);
#                 end
#                 zs=zeros(1,nZeros);
#                 for i = 1:nZeros,
#                     line=fgetl(fid);
#                     val=sscanf(line,'%f %f');
#                     zs(i)=val(1)+j*val(2);
#                 end
#                 ts_verbose(0,'PAZ2 line: %d Poles, %d Zeros, "%s"',nPoles,
#                            nZeros, descript);
#                 pazDone=1;
#             case 'DIG2',    % Digitizer response (only one) [Same for GSE 2.0
#                             %and 2.1]
#                 if digDone, error('Multiple DIG2 lines'); end
#                 %sNum=str2double(line(6:7));    % Stage sequence number
#                 digGain=str2double(line(9:23));
#                 inSPS=  str2double(line(25:35));
#                 descrip=strtrim(line(37:end));
#                 ts_verbose(0,'DIG2 line: input samprate = %g, digitizer type"
#                            " = ''%s''', inSPS,descrip);
#                 digDone=1;
#             otherwise
#                 ts_warning('Unknown GSE2 line: ''%s''',line(1:4));
#         end
#         line=fgetl(fid);
#         if ~ischar(line), break, end
#     end
#
#     % Calculate and verify gain
#     % We want gain in counts/nm, digGain is counts/V, sensNormConst is V/nm
#     % Check using refGain (nm/count) at refFreq
#     gain = digGain*sFactor;
#     tester=ts_response('PZ','nm','counts',gain,ps,zs);
#     testGain=abs(tester.response(1/calPer));
#     if testGain*calib >1.01 || testGain*calib < 0.99;
#         ts_warning('Prod of stage gains (%g nm/c) dis w/ calib gain'
#                    ' (%g nm/c) at ref freq (%g Hz)', ...
#             1/testGain,calib,1/calPer);
#     end
#
#     % Fill the GSEprops structure
#     %onDate
#     %onTime
#     GSEprops.startDate  = ts_datetime([onDate ' ' onTime]);
#     GSEprops.refFreq    = 1/calPer;
#     GSEprops.station    = station;
#     GSEprops.comp       = comp;
#     GSEprops.inSPS      = inSPS;
#     GSEprops.outSPS     = outSPS;
#     GSEprops.sensorType = instType;
#     GSEprops.digType = descrip;
#     GSEprops.digGain = digGain;
#
# end
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function t = cleanUpOld(t,inPZ,nPZ,iPZ)
#     if inPZ==0; return; end
#     while iPZ < nPZ,
#         if inPZ==1,
#             t.poles=[t.poles 0];
#         elseif inPZ==2,
#             t.privateZeros=[t.privateZeros 0];
#         else
#             error('invalid inPZ value: %d',inPZ);
#         end
#         iPZ = iPZ+1;
#     end
# end
#
