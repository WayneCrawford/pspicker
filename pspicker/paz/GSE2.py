import copy

import numpy as np
from obspy.core.inventory.response import PolesZerosResponseStage

from .PAZ import _get_paz


def write_GSE2(fname, sta_code, ch, debug=False):
    """
    Write channel to the GSE2 format used by SEISAN

    GSE2 specification is found in GSE_2.0.pdf
    :param fname: name of file to write to
    :param sta_code: station code
    :param ch: obspy Channel objet
    """
    if debug:
        print(f'writing to {fname}\n')
        ch.response.plot(0.1, 'VEL')
    with open(fname, 'w') as f:
        f.write(_CAL2_line(ch, sta_code))
        f.write(_PAZ2_lines(ch))
        f.write(_DIG2_line(ch))


def _CAL2_line(ch, sta_code):
    """
    Returns Calibration Identification Line (CAL2)

    From Table 13, page 69
    Gives a rapid idea of sensitivity in passband
    """
    id = 'CAL2'.ljust(4)[:4]
    auxid = ''
    inst_type = ch.sensor.model
    # system gain (count/m) at ref freq (calfreq)
    calib, calfreq = _get_gain_disp(ch.response)
    print('CAL2 gain = {:.3g}'.format(abs(calib)))
    sens = 1.e9 / abs(calib)  # nm/c
    samprat = ch.sample_rate
    cal_line = '{} {} {} {} {} {:10.2e} {:7.3f} {:10.5f} {} {}\n'.format(
        id.ljust(4)[:4], sta_code.ljust(5)[:5], ch.code.ljust(3)[:3],
        auxid.ljust(4)[:4], inst_type.ljust(6)[:6],
        sens, 1 / calfreq, samprat,
        ch.start_date.strftime('%Y/%m/%d %H:%M'),
        ch.end_date.strftime('%Y/%m/%d %H:%M'))
    return cal_line


def _PAZ2_lines(ch):
    """
    Returns Pole-zero Response Section (PAZ2)

    From Table 14, page 70
    """
    paz = _get_paz(ch.response)
    # id = 'PAZ2'
    snum = paz.stage_sequence_number
    print(f'Last PAZ Stage number = {snum:d}')
    for st in ch.response.response_stages[0: snum]:
        print(st)
    # gain is in V/m
    gain, freq = _get_gain_disp(ch.response, 1, snum)
    # change poles/zeros to V/m also
    poles, zeros = _vel_to_disp(paz.poles, paz.zeros)
    nfact = _pz_norm_fact(poles, zeros, freq)
    gain_norm = gain * nfact
    print('PAZ gain = {abs(gain):.3g}')
    print('pz_norm_fact = {nfact:.3g}')
    norm_const = abs(gain_norm)/1e9
    deci = ''
    corr = ''
    npole = len(poles)
    nzero = len(zeros)
    descrip = paz.pz_transfer_function_type
    paz_header = f'PAZ2 {snum:2d} V {norm_const:15.8e} {deci.ljust(4)[:4]} '
    paz_header += f'{corr.ljust(8)[:8]} {npole:3d} {nzero:3d} '
    paz_header += f'{descrip.ljust(25)[:25]}\n'
    p_lines = ''
    for val in poles:
        p_lines += f' {val.real:15.8e} {val.imag:15.8e}\n'
    z_lines = ''
    for val in zeros:
        z_lines += f' {val.real:15.8e} {val.imag:15.8e}\n'
    return paz_header + p_lines + z_lines


def _DIG2_line(ch):
    """
    Returns Digitizer Response Section (DIG2)

    Assumes digitizer is right after the PZ stage(s)
    From Table 17, page 72
    """
    stages = [copy.deepcopy(stage) for stage in ch.response.response_stages
              if not isinstance(stage, PolesZerosResponseStage)]
    digi_stage = stages[0]
    id = 'DIG2'
    snum = digi_stage.stage_sequence_number
    sensitivity = digi_stage.stage_gain   # counts/input unit
    print(f'DIG2 gain = {sensitivity:.3g}')
    samprat = digi_stage.decimation_input_sample_rate
    descrip = digi_stage.description.replace('DIGITIZER', '').lstrip(' -')
    digi_line = f'{id.ljust(4)[:4]} {snum:2d} {sensitivity:15.8e} '
    digi_line += f'{samprat:11.5f} {descrip.ljust(25)[:25]}\n'
    # digi_line = id.ljust(4)[:4]  + ' ' + \
    #             '{:2d}'.format(snum)  + ' ' + \
    #             '{:15.8e}'.format(sensitivity)  + ' ' + \
    #             '{:11.5f}'.format(samprat)  + ' ' + \
    #             descrip.ljust(25)[:25]  + '\n'
    return digi_line


def _vel_to_disp(poles, zeros):
    """
    Convert velocity pole-zeros to displacement pole-zeros

    Removes a zero at zero frequency, if there is one, otherwise adds
    a poles at zero frequency
    """
    new_zeros = zeros.copy()
    new_poles = poles.copy()
    new_zeros.append(0. + 0.j)
    return new_poles, new_zeros


def _get_gain_disp(response, start_stage=None, end_stage=None):
    """
    Return stage(s) displacement gain
    """
    freq = response.instrument_sensitivity.frequency
    # resp is in outputs/inputs (full instrument = counts/m)
    gain = response.get_evalresp_response_for_frequencies(
        [freq], 'DISP', start_stage, end_stage)
    return gain[0], freq


def _pz_norm_fact(poles, zeros, ref_freq):
    """
    Returns pole-zero normalisation factor

    poles and zeros expressed in radians
    ref_freq in Hz
    """
    norm_fact = 1.
    for p in poles:
        norm_fact *= (2*np.pi*ref_freq*1j - p)
        print(abs(norm_fact))
    for z in zeros:
        norm_fact /= (2*np.pi*ref_freq*1j - z)
        print(abs(norm_fact))
    norm_fact = abs(norm_fact)
    return norm_fact
