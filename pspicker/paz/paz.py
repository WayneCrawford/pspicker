"""
Pole-Zero representation of an instrument

Follows SEED Pole-Zero protocol for Analog Stages, radian format
(Seed Manual v2.4, Appendix C)

G(f) = S_d * A_0 * prod(s - z_n) / prod(s - p_n) = S_d * A_0 * Hp(s)
where s=j*2*pi*f

The corresponding attributes are:
    S_d: gain factor
    A_0: norm_factor
    G(f): gain at frequency f
    ref_gain = G(f0) = gain at the ref_freq. If there is a passband (a place
                       where the response is flat, ref_freq should be chosen
                       in this band)

The user can specify gain or ref_gain, not both
If norm_factor is not specified, it is automatically calculated so that
    A_0 * abs(Hp(s)) = 1 at f0
as per the SEED recommendation

The user can also change the input_units after creating the class, in which
case the poles, zeros and gain will  be modified if the starting and ending
units are in self.known_units().
"""
import numpy as np
import json
import copy
import warnings

from obspy.core.inventory import read_inventory
from obspy.core.inventory.response import PolesZerosResponseStage
from matplotlib import pyplot as plt


class PAZ():
    """
    Object-based implementation of Pole-Zero representation

    For representing a full instrument for simulations, etc
    """
    converter = {"m":       [1.,    2],
                 "m/s":     [1.,    1],
                 "m/s^2":   [1.,    0],
                 "m/s**2":  [1.,    0],
                 "mm":      [1.e-3, 2],
                 "mm/s":    [1.e-3, 1],
                 "mm/s^2":  [1.e-3, 0],
                 "mm/s**2": [1.e-3, 0],
                 "nm":      [1.e-9, 2],
                 "nm/s":    [1.e-9, 1],
                 "nm/s^2":  [1.e-9, 0],
                 "nm/s**2": [1.e-9, 0],
                 "gal":     [1.e-2, 0],   # 1 cm/s^2
                 "mgal":    [1.e-5, 0],
                 "ugal":    [1.e-8, 0]}

    def __init__(self, gain, poles=[], zeros=[], ref_freq=1.,
                 input_units='m/s', output_units='counts',
                 norm_factor=None):
        """
        Initialize using pole-zero gain (NOT passband gain)

        Pole and zero Frequency units are radians, ref_freq is Hertz

        :param gain: pole-zero gain (NOT passband gain)
        :param poles: poles (rad/s)
        :param zeros: zeros (rad/s)
        :param ref_freq: reference frequency (Hz) for passband gain
        :param input_units: units (usually physical) at the sensor entry
        :param output_units: units as stored
        :param norm_fact: A0 normalization factor
        """
        self.gain = gain
        self.poles = np.array(poles)
        self.zeros = np.array(zeros)
        self.ref_freq = ref_freq
        self._input_units = input_units
        self.output_units = output_units
        if norm_factor is None:
            self._norm_factor = self.calc_norm_factor()
        else:
            self._norm_factor = norm_factor

    @classmethod
    def from_refgain(cls, ref_gain, poles=[], zeros=[], ref_freq=1.,
                     input_units='m/s', output_units='counts',
                     norm_factor=None):
        """
        Initialize using gain at ref_freq

        Pole and zero Frequency units are radians, ref_freq is Hertz

        :param ref_gain: gain at ref_freq (should be w/in instrument passband)
        :param poles: poles (rad/s)
        :param zeros: zeros (rad/s)
        :param ref_freq: reference frequency (Hz) for passband gain
        :param input_units: units (usually physical) at the sensor entry
        :param output_units: units as stored
        :param norm_fact: A0 normalization factor
        """
        cls = PAZ(1, poles=poles, zeros=zeros, ref_freq=ref_freq,
                  input_units=input_units,
                  output_units=output_units,
                  norm_factor=norm_factor)
        cls.gain = ref_gain * cls.calc_norm_factor() / cls.norm_factor
        return cls

    @property
    def ref_gain(self):
        # return self.gain / self.norm_factor
        return self.gain * self.norm_factor / self.calc_norm_factor()

    @property
    def norm_factor(self):
        """ A0 normalization factor """
        return self._norm_factor

    @norm_factor.setter
    def norm_factor(self, value):
        """gain must change if norm factor changed"""
        self.gain *= self._norm_factor / value
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
        try:
            self._unit_conversion(value)
            self._input_units = value
        except Exception:
            print('Could not convert input_units to "{}", keeping "{}"'
                  .format(value, self._input_units))

    def _hp(self, freqs):
        """
        Return Hp(f) = prod(s - z_n) / prod(s - p_n)

        where s=j*2*pi*f
        """
        s = 2 * np.pi * np.array(freqs) * 1j
        hp = np.ones(s.shape, np.complex128)
        for p in self.poles:
            hp /= (s - p)
        for z in self.zeros:
            hp *= (s - z)
        return hp

    def __eq__(self, obj, places=2, verbose=False):
        """
        Return true if objects are equal to within places decimal places

        :param places: decimal place precision used for gain
        """
        err_msgs = []
        diff = abs((self.gain*self.norm_factor)
                   / (obj.gain*obj.norm_factor) - 1) * 10**places
        if diff > 1:
            err_msgs.append('different gain*A0 ()')
        if not len(self.zeros) == len(obj.zeros):
            err_msgs.append('n_zeros different')
        if not len(self.poles) == len(obj.poles):
            err_msgs.append('n_poles different')
        for z1, z2 in zip(self.zeros, obj.zeros):
            if not z1 == z2:
                err_msgs.append('zero different')
        for p1, p2 in zip(self.poles, obj.poles):
            if not p1 == p2:
                err_msgs.append('pole different')
        if not self.input_units.upper() == obj.input_units.upper():
            err_msgs.append('different input units')
        if not self.output_units.upper() == obj.output_units.upper():
            err_msgs.append('different output units')
        if not self.ref_freq == obj.ref_freq:
            err_msgs.append('different ref_freqs')
        if len(err_msgs) > 0:
            if verbose:
                print('PAZs NOT EQUAL:' + ','.join(err_msgs))
            return False
        return True

    def __str__(self):
        s = 'PAZ: '
        s += f'gain={self.gain:.4g}, '
        s += f'A0={self.norm_factor:.4g}, '
        s += f'poles={self.poles}, '
        s += f'zeros={self.zeros}, '
        s += f'(ref_gain={self.ref_gain:.3g} '
        s += f'{self.output_units}/{self.input_units} '
        s += f'at {self.ref_freq:.3g}Hz)'
        return s

    def copy(self):
        return copy.deepcopy(self)

    def get_response(self, freqs):
        """ return response for a given range of frequencies"""
        resp = self.gain * self.norm_factor * self._hp(freqs)
        # s = 2 * np.pi * freqs * 1j
        # resp = self.gain * self.norm_factor * np.ones(s.shape, np.complex128)
        # for p in self.poles:
        #     resp /= (s - p)
        # for z in self.zeros:
        #     resp *= (s - z)
        return resp

        s = 2 * np.pi * freqs * 1j
        hp = np.ones(s.shape, np.complex128)
        for p in self.poles:
            hp /= (s - p)
        for z in self.zeros:
            hp *= (s - z)
        return hp

    def plot(self, min_freq, output='VEL', show=True, label='', axes=None,
             sampling_rate=None, sym='b-'):
        """
        Plot PAZ object's Response

        :type min_freq: float
        :param min_freq: Lowest frequency to plot.
        :type output: str
        :param output: Output units. One of:

                ``"DISP"``
                    displacement
                ``"VEL"``
                    velocity
                ``"ACC"``
                    acceleration
        :type axes: list of 2 :class:`matplotlib.axes.Axes`
        :param axes: List/tuple of two axes instances to plot the
            amplitude/phase spectrum into. If not specified, a new figure is
            opened.
        :type label: str
        :param label: Label string for legend.
        :type sampling_rate: float
        :param sampling_rate: Manually specify sampling rate of time series.
            If not given it is attempted to determine it from the information
            in the individual response stages.  Does not influence the spectra
            calculation, if it is not known, just provide the highest frequency
            that should be plotted times two.
        :type show: bool
        :param show: Whether to show the figure after plotting or not. Can be
            used to do further customization of the plot before showing it.
        :type sym: str
        :param sym: symbol to use on plot
        """
        paz = self.copy()
        if output == 'DISP':
            paz.input_units = 'm'
        elif output == 'VEL':
            paz.input_units = 'm/s'
        elif output == 'ACC':
            paz.input_units = 'm/s^2'
        else:
            print(f'Unknown output type: {output}, using VEL')
            paz.input_units = 'm/s'
        if sampling_rate is None:
            if min_freq < 1:
                sampling_rate = 100.
            else:
                sampling_rate = min_freq * 100.
        else:
            assert sampling_rate > 2*min_freq
        if axes is not None:
            axs = axes
            fig = axs[0].figure
            npts = 50
        else:
            fig, axs = plt.subplots(2, 1, sharex=True)
            npts = 100
        f = np.power(10., np.linspace(np.log10(min_freq),
                                      np.log10(sampling_rate/2), npts))
        resp = paz.get_response(f)
        axs[0].loglog(f, np.absolute(resp), sym, label=label)
        axs[0].set_ylabel('Amplitude ({}/{})'.format(self.output_units.lower(),
                                                     self.input_units.lower()))
        axs[0].grid(True)
        axs[0].legend()
        axs[1].semilogx(f, np.angle(resp), sym)
        axs[1].grid(True)
        axs[1].set_ylabel('Phase')
        axs[1].yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
        axs[1].yaxis.set_major_formatter(plt.FuncFormatter(
            multiple_formatter()))
        axs[1].set_xlabel('Frequency (Hz)')
        if show is True:
            plt.show()

    def calc_norm_factor(self):
        """
        Calculate A0 normalization factor

        Returns A0 such that A0 * prod(s-z_n)/prod(s-p_n) = 1 at ref_freq
        (SEED recommendation)
        """
        hp = self._hp([self.ref_freq])
        val = 1./abs(hp[0])
        # s = 2 * np.pi * self.ref_freq * 1j
        # val = 1
        # for p in self.poles:
        #     val *= (s - p)
        # for z in self.zeros:
        #     val /= (s - z)
        # val = abs(val)
        return val

    def to_obspy(self):
        """
        Return obspy paz dict

        obspy is not clear about what is in the "gain" and "sensitivity"
        keys.  I tried self.norm_factor and self.gain, respectively, but
        this makes their product depend on f_ref.  Switched to:
            self.norm_factor, self.ref_gain
        """
        return {'poles': list(self.poles), 'zeros': list(self.zeros),
                'f_ref': self.ref_freq, 'gain': self.norm_factor,
                'sensitivity': self.ref_gain}

    def _unit_conversion(self, new_input_units):
        """
        Convert from self.input_units to new_input_units
        """
        igain, inz = self._ref_meter(self.input_units)
        ogain, onz = self._ref_meter(new_input_units)

        # print(f'{self.input_units=}, {igain=}, {inz=}')
        # print(f'{new_input_units=}, {ogain=}, {onz=}')

        if igain is None or ogain is None:
            raise ValueError('Could not convert')
            # warnings.warn('Did not convert')
            # return

        self.gain *= ogain / igain
        # if self._norm_factor is not None:
        #     old_calc_norm_factor = self.calc_norm_factor()
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
        # if self._norm_factor is not None:
        #     # Since norm factor cannot change, change gain
        #     self.gain *= old_calc_norm_factor / self.calc_norm_factor()

    @classmethod
    def from_obspy_PoleZeroResponseStage(cls, stage):
        """
        Return PAZ correpsonding to obspy PoleZeroResponseStage
        """
        assert isinstance(stage, PolesZerosResponseStage)
        assert stage.stage_gain_frequency == stage.normalization_frequency

        # sampling_rate = None
        # if stage.decimation_input_sample_rate is not None:
        #     sampling_rate = stage.decimation_input_sample_rate
        #     if stage.decimation_factor is not None:
        #         sampling_rate /= stage.decimation_factor

        cls = PAZ(stage.stage_gain, poles=stage.poles, zeros=stage.zeros,
                  ref_freq=stage.stage_gain_frequency,
                  input_units=stage.input_units,
                  output_units=stage.output_units,  
                  norm_factor=stage.normalization_factor)
        return cls

    @classmethod
    def from_obspy_response(cls, resp):
        """
        Return PAZ fitting an obspy response object

        Combines all poles and zeros, multiplies all gains
        Does not fit frequency responses of other stage types
        """
        # print(f'{resp=}')
        obspy_paz = _get_pzrs_all(resp)
        # print(f'{obspy_paz=}')
        cls = PAZ.from_obspy_PoleZeroResponseStage(obspy_paz)
        return cls

    def write_json_pz(self, filename):
        """
        write JSON paz format
        """
        if not self.output_units.lower() == 'counts':
            print("Can't write to json_pz: output_units are "
                  f"{self.output_units}, not counts")
            return
        out_paz = self.copy()
        out_paz.input_units = 'm/s'
        # out_paz.ref_freq=.1
        outdict = {'poles': [[x.real, x.imag] for x in out_paz.poles],
                   'zeros': [[x.real, x.imag] for x in out_paz.zeros],
                   'sensitivity': 1/out_paz.ref_gain,
                   'input_units': out_paz.input_units,
                   'norm_const': out_paz.norm_factor,
                   'f_ref': out_paz.ref_freq}
        # if out_paz.sampling_rate is not None:
        #     outdict['sampling_rate'] = out_paz.sampling_rate
        with open(filename, 'w') as f:
            try:
                json.dump(outdict, f, indent=2)
            except Exception as e:
                raise ValueError(f'Could not write to {filename}: {e}')

    @classmethod
    def read_json_pz(cls, filename):
        """
        read JSON paz format
        """
        with open(filename) as f:
            try:
                paz = json.load(f)
            except Exception as e:
                raise ValueError(f'Could not load {filename}: {e}')
        poles, zeros = [], []
        for p in paz['poles']:
            poles.append(float(p[0]) + float(p[1])*1.j)
        for z in paz['zeros']:
            zeros.append(float(z[0]) + float(z[1])*1.j)
        val = cls.from_refgain(1./float(paz['sensitivity']),
                               poles=poles, zeros=zeros,
                               ref_freq=float(paz['f_ref']),
                               input_units=paz.get('input_units', 'm/s'),
                               output_units='counts',
                               norm_factor=paz.get('norm_const', None))
        return val

    @classmethod
    def read_stationxml(cls, filename, channel, station='*'):
        """
        read JSON paz format
        """
        inv = read_inventory(filename, 'STATIONXML')
        inv = inv.select(channel=channel, station=station)
        nscs = [{'net': n.code, 'sta': s.code, 'chan': c}
                for n in inv for s in n for c in s]
        if len(nscs) > 1:
            nsc_strs = [f'{x["net"]}.{x["sta"]}.{x["chan"].code}'
                        for x in nscs]
            print('Multiple channels were selected, using {nsc_strs[0]}')
            print('Other channels were: ' + ','.join(nsc_strs[1:]))
        elif len(nscs) == 0:
            raise NameError(f'Did not find {station}, {channel} in {filename}')
        resp = nscs[0]["chan"].response
        return cls.from_obspy_response(resp)

    @classmethod
    def read_sac_pz(cls, filename):
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
        # constant is norm_factor * gain, so if I set one to constant, I have
        # to set the other to 1.
        val = cls(constant,
                  poles=poles,
                  zeros=zeros,
                  norm_factor=1,
                  input_units='m',
                  output_units='counts')
        return val

    @classmethod
    def read_gse_response(cls, filename):
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
        val = cls.from_refgain(1./info['sensitivity'],
                               poles=info['poles'],
                               zeros=info['zeros'],
                               ref_freq=info['f_ref'],
                               input_units='nm',
                               output_units='counts',
                               norm_factor=info.get('norm_const', None))
        return val

    @classmethod
    def read_baillard_pz(cls, filename):
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
        val = cls.from_refgain(1/info['sensitivity'],
                               poles=info['poles'],
                               zeros=info['zeros'],
                               ref_freq=info['f_ref'],
                               input_units='nm',
                               output_units='counts',
                               norm_factor=info['A0'])
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
            warnings.warn('Unknown conversion unit: {}, known values are "{}"'
                          .format(units, ', '.join(self.converter.keys())))
            return None, None


def _get_pzrs(resp, verbose=False):
    """
    Get PolesZerosResponseStages

    Like obspy get_paz(), except that. if there is more than one PZ stage,
    combines them into one

    Does not handle decimation, resource_id, resource_id2
    Returns LAST stage's stage_sequence number
    :param resp: obspy Response object
    :returns: obspy PoleZeroResponse stage (stage_sequence_number = LAST
        paz stage)
    """

    pz_stages = [copy.deepcopy(stage) for stage in resp.response_stages
                 if isinstance(stage, PolesZerosResponseStage)]
    if len(pz_stages) == 0:
        print("No PolesZerosResponseStage found.")
        pz_out = PolesZerosResponseStage(normalization_factor=1.,
                                         stage_gain=1.,
                                         poles=[],
                                         zeros=[])
    elif len(pz_stages) == 1:
        pz_out = pz_stages[0]
    else:
        if verbose:
            print(f"Combining {len(pz_out):d} PolesZerosResponseStages")
        pz_out = pz_stages[0]
        for pz in pz_stages[1:]:
            if _is_compatible(pz_out, pz):
                if pz.description is not None:
                    pz_out.description += ' + ' + pz.description
                if not pz_out.name:
                    if pz.name:
                        pz_out.name = ' + ' + pz.name
                elif pz.name:
                    pz_out.name += ' + ' + pz.name
                pz_out.normalization_factor *= pz.normalization_factor
                pz_out.output_units = pz.output_units
                pz_out.output_units_description = pz.output_units_description
                pz_out.stage_gain *= pz.stage_gain
                pz_out.poles.extend(pz.poles)
                pz_out.zeros.extend(pz.zeros)
                pz_out.stage_sequence_number = pz.stage_sequence_number
            else:
                print('Incompatible stages... quitting')
                return None
    return pz_out


def _get_pzrs_all(resp):
    """
    Get PoleZeroResponseStage covering ALL stages

    :param resp: obspy Response object
    :returns: obspy PoleZeroResponse stage (stage_sequence_number = LAST
        PoleZeroResponseStage)
    """
    pz_out = _get_pzrs(resp)
    for stage in resp.response_stages[pz_out.stage_sequence_number:]:
        if _is_compatible(pz_out, stage, next_is_PAZ=False):
            pz_out.output_units = stage.output_units
            pz_out.output_units_description = stage.output_units_description
            pz_out.stage_gain *= stage.stage_gain
            pz_out.stage_sequence_number = stage.stage_sequence_number
            pz_out.decimation_input_sample_rate =\
                stage.decimation_input_sample_rate
            pz_out.decimation_factor = stage.decimation_factor
        else:
            print('Incompatible stages... quitting')
            return None
    if pz_out.description is None:
        pz_out.description = 'PAZ stages'
    pz_out.description += ' + non-PAZ stages'
    if pz_out.name:
        pz_out.name += ' + non-PAZ stages'
    return pz_out


def _is_compatible(first_stage, next_stage, next_is_PAZ=True):
    """Can one PAZ stage be safely added to another?"""
    if not first_stage.output_units == next_stage.input_units:
        msg = ("First stage out units not equal to next stage in units. ")
        warnings.warn(msg)
        return False
    if next_is_PAZ:
        if not (first_stage.pz_transfer_function_type
                == next_stage.pz_transfer_function_type):
            msg = ("Stage transfer function types not equal")
            warnings.warn(msg)
            return False
        if not (first_stage.stage_gain_frequency
                == next_stage.stage_gain_frequency):
            print("Stage gain frequencies not equal ({:g} and {:g})".format(
                first_stage.stage_gain_frequency,
                next_stage.stage_gain_frequency))
            return False
        if not (first_stage.normalization_frequency
                == next_stage.normalization_frequency):
            msg = ("Normalization frequencies not equal")
            warnings.warn(msg)
            return False
        if first_stage.decimation_factor or next_stage.decimation_factor:
            msg = ("Can't handle decimation (yet)")
            warnings.warn(msg)
            return False
    return True


def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
    def gcd(a, b):
        while b:
            a, b = b, a % b
        return a

    def _multiple_formatter(x, pos):
        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num, den)
        (num, den) = (int(num/com), int(den/com))
        if den == 1:
            if num == 0:
                return r'$0$'
            if num == 1:
                return r'$%s$' % latex
            elif num == -1:
                return r'$-%s$' % latex
            else:
                return r'$%s%s$' % (num, latex)
        else:
            if num == 1:
                return r'$\frac{%s}{%s}$' % (latex, den)
            elif num == -1:
                return r'$\frac{-%s}{%s}$' % (latex, den)
            else:
                return r'$\frac{%s%s}{%s}$' % (num, latex, den)
    return _multiple_formatter
