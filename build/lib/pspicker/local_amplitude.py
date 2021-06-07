# import os.path
#
# import json

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from scipy.signal import find_peaks
from obspy.core.event.magnitude import Amplitude
from obspy.core.event.origin import Pick
from obspy.signal.invsim import estimate_wood_anderson_amplitude

from .logger import log
from .paz import PAZ


class Amp():
    """Mini class to hold and print amplitude values"""
    def __init__(self, value=0, period=0, time=None, channel=None):
        self.value = value
        self.period = period
        self.time = time
        self.channel = channel

    def __str__(self):
        return "Amp: value={:.3g}, period={:.3g}, time={}, chan={}".format(
            self.value, self.period, self.time, self.channel)


class LocalAmplitude():
    """
    Calculate local earthquake signal amplitude at a station
    """
    def __init__(self, traces, picks, response_file, response_file_type,
                 P_pick_window=[-5, 30], S_pick_window=[-20, 10]):
        """
        :param traces: all traces associated with this station/event
        :param picks: P- and S- Picks
        :type picks: list of obspy Pick
        :param response_file: file containing instrument response
        :param response_file_type:
        :param P_pick_window: if there is only a P pick, window before and
            after to look for maximum
        :param P_pick_window: if there is only an s pick, window before and
            after to look for maximum
        """
        assert len(picks) > 0, "no picks!"
        for pick in picks:
            assert isinstance(pick, Pick), 'pick is not an obspy Pick'
        self.pick_P = self._get_pick(picks, "P")
        self.pick_S = self._get_pick(picks, "S")
        self.P_pick_window = P_pick_window
        self.S_pick_window = S_pick_window
        self.win_start, self.win_end, self.ref_pick, self.station =\
            self._set_ampl_window()
        # Should take advantage of obspy to read all standard formats
        self.paz = get_response(response_file, response_file_type)
        self.paz.input_units = 'nm'
        self.traces = traces

    def __str__(self):
        """
        """
        s = "LocalAmplitude: "
        if self.pick_P:
            s += f"P_time = {self.pick_P.time}, "
        if self.pick_S:
            s += f"S_time = {self.pick_S.time}, "
        if self.win_start:
            s += f"win_start = {self.win_start}, "
        if self.win_end:
            s += f"win_end = {self.win_end}, "
        if self.traces:
            s += "traces = ({}), ".format(
                ",".join([x.stats.channel for x in self.traces]))
        s += "paz = {}, ".format(self.paz)
        s += "P_pick_window = {}, ".format(self.P_pick_window)
        s += "S_pick_window = {}".format(self.S_pick_window)
        return s

    def get_iaml(self, plot=False, method='wood_calc'):
        """
        get IAML amplitude and associated pick
        From IAPEI CoSOI 2013 Working Group recommentations:
        "ML = log10(A) + 1.11log10(R) + 0.00189*R - 2.09
        where
        A = maximum *trace* amplitude in *nm* that is measured on output from
            a *horizontal-component* instrument that is filtered so that the
            response of the seismograph/filter system replicates that of a
            *Wood-Anderson standard seismograph* but with a
            static magnification of 1
        R = *hypocentral distance in km*, typically less than 1000 km

        the constant -2.09 is based on an experimentally determined static
        magnification of the Wood-Anderson of 2080, rather than the theoretical
        magnification of 2800.

        The amplitudes used in the magnitude formulas ... are in most
        circumstances to be measured as one-half the maximum deflection of the
        seismogram trace, peak-to-adjacent- trough or trough-to-adjacent-peak,
        where peak and trough are separated by one crossing of the zero-line:
        the measurement is sometimes described as “one-half peak-to-peak
        amplitude.” None of the magnitude formulas presented in this article
        are intended to be used with the full peak-to-trough deflection as
        the amplitude.

        WA displacement zeros = [(0+0j), (0+0j)]
        WA displacement poles = [(-5.49779 - 5.60886j), (-5.49779 + 5.60886j)]
        WA seismometer free period = 0.8 s
        WA damping constant = 0.7
        A0 = 0.97866 @ 4 Hz
        :param plot: plot the result
        :param method:
            "wood_calc": calculate amplitude on the signal transformed to
                         simulate a unity gain Wood-Anderson filter
            "wood_est": calculate amplitude on the original signal and use
                        obspy.signal.invsim.estimate_wood_anderson_amplitude()
            "raw_disp": calculate amplitude on the signal transformed to
                        displacement in nm
        :returns: amplitude in meters (type=IAML), associated Pick
        :rtype: ~class obspy.core.event.magnitude.Amplitude,
                ~class obspy.core.event.origin.Pick
        """
        if self.ref_pick is None:
            log('No ref_pick!', 'debug')
            return None
        assert method in ['wood_calc', 'wood_est', 'raw_disp']

        signal = self.traces.copy()
        # do all work in nm
        paz_simulate, paz_simulate_obspy = None, None
        paz_remove = self.paz.copy()
        if method == 'wood_est':
            paz_remove.input_units = 'm/s'
            plot_units = 'Original (counts)'
        else:
            if method == 'wood_calc':
                paz_simulate = PAZ.from_refgain(
                    1,
                    poles=[(-5.49779 - 5.60886j), (-5.49779 + 5.60886j)],
                    zeros=[(0+0j), (0+0j)],
                    ref_freq=4.0,
                    input_units='nm', output_units='counts')
                plot_units = 'Wood-And (nm)'
                paz_simulate_obspy = paz_simulate.to_obspy()
            elif method == 'raw_disp':
                paz_simulate = None
                paz_remove.input_units = 'nm'
                plot_units = 'Disp (nm)'
            for tr in signal:
                tr.simulate(paz_remove=paz_remove.to_obspy(),
                            paz_simulate=paz_simulate_obspy,
                            water_level=60.0)
        amp = pk2pk(signal, self.win_start, self.win_end)
        if amp is None:
            return None, None
        if method == 'wood_est':
            amp.value = estimate_wood_anderson_amplitude(paz_remove.to_obspy(),
                                                         amp.value, amp.period)
            # estimate...() gives displacement amplitude ON a Wood-Anderson
            # seismometer, divide by 2080 then multiply by 1e6 to get ground
            # motion in nm?
            amp.value /= 2080./1e6
        if plot:
            self.plot(signal, amp, plot_units)
        waveform_id = self.ref_pick.waveform_id
        waveform_id.channel_code = amp.channel
        pick = Pick(time=amp.time,
                    waveform_id=waveform_id,
                    method_id='ps_picker',
                    phase_hint='IAML',
                    evaluation_mode='automatic',
                    evaluation_status='preliminary')
        obspy_amp = Amplitude(generic_amplitude=amp.value/(2. * 1.e9),
                              type='IAML',
                              unit='m',
                              period=amp.period,
                              magnitude_hint='ML',
                              category='period',
                              pick_id=pick.resource_id,
                              waveform_id=pick.waveform_id)
        return obspy_amp, pick

    def _set_ampl_window(self):
        if self.pick_S is not None and self.pick_P is not None:
            pick = self.pick_S
            win_start = self.pick_P.time - 1
            win_end = self.pick_S.time + 2*(self.pick_S.time
                                            - self.pick_P.time)
        elif self.pick_S is not None:
            pick = self.pick_S
            win_start = self.pick_S.time + self.S_pick_window[0]
            win_end = self.pick_S.time + self.S_pick_window[1]
        elif self.pick_P is not None:
            pick = self.pick_P
            win_start = self.pick_P.time + self.P_pick_window[0]
            win_end = self.pick_P.time + self.P_pick_window[1]
        else:
            return None, None, None, None
        return win_start, win_end, pick, pick.waveform_id.station_code

    @staticmethod
    def _get_pick(picks, phase):
        pick = [x for x in picks if x.phase_hint[0] == phase]
        if len(pick) == 0:
            return None
        else:
            if not len(pick) == 1:
                log('{:d} picks have phase[0] == "{}": {}, returning first one'
                    .format(len(pick), pick[0].phase_hint[0], pick), 'error')
            return pick[0]

    def plot(self, transformed, Amp, trans_units):
        """
        Plot the amplitude pick

        Top plot: all traces, colored by type
        Middle trace: Wood-Anderson traces, colored by type, with ampl pick
        Bottom trace: zoom on amplitude region, with max pick for max trace
        :param transformed: transformed traces
        :param Amp: pk2pk output
        :param trans_units: y-axis label for the transformed traces
        """
        fig, axs = plt.subplots(3, 1, num=f'LocalAmplitude {self.station}')
        traces = self.traces.slice(self.win_start, self.win_end)
        transcut = transformed.slice(self.win_start, self.win_end)
        imin = transcut[0].times('utcdatetime').searchsorted(Amp.time)
        imax = transcut[0].times('utcdatetime')\
                          .searchsorted(Amp.time + Amp.period) - 1
        zoom_start = Amp.time - 2*Amp.period
        zoom_end = Amp.time + 4*Amp.period

        # Top axis: plot all traces
        axs[0].set_xlim(self.win_start.matplotlib_date,
                        self.win_end.matplotlib_date)
        for tr in traces:
            color, label = self._choose_color(tr)
            axs[0].plot_date(tr.times(type="matplotlib"), tr.data, color + '-',
                             label=label)
        axs[0].add_patch(Rectangle((zoom_start.matplotlib_date,
                                    axs[0].viewLim.y0),
                                   width=(zoom_end - zoom_start)/86400,
                                   height=axs[0].viewLim.y1-axs[0].viewLim.y0,
                                   fc='blue', alpha=0.5, zorder=-1))
        if self.pick_P:
            axs[0].axvline(self.pick_P.time.matplotlib_date, color='b')
        if self.pick_S:
            axs[0].axvline(self.pick_S.time.matplotlib_date, color='r')
        axs[0].set_ylabel(self.station)
        axs[0].legend(loc='center left')

        # Middle axis: plot transformed traces
        axs[1].clear()
        axs[1].set_xlim(self.win_start.matplotlib_date,
                        self.win_end.matplotlib_date)
        for tr in transcut:
            color, label = self._choose_color(tr)
            axs[1].plot_date(tr.times(type="matplotlib"), tr.data, color + '-',
                             label=label)
        axs[1].add_patch(Rectangle((zoom_start.matplotlib_date,
                                    axs[1].viewLim.y0),
                                   width=(zoom_end - zoom_start)/86400,
                                   height=axs[1].viewLim.y1-axs[1].viewLim.y0,
                                   ec='blue', fc='lightblue', alpha=0.3,
                                   zorder=-1))
        axs[1].set_ylabel(trans_units)
        axs[1].legend(loc='center left')

        axs[2].set_xlim(zoom_start.matplotlib_date, zoom_end.matplotlib_date)
        tr = transcut.select(channel=Amp.channel)[0]
        color, label = self._choose_color(tr)
        axs[2].plot_date(tr.times(type="matplotlib"), tr.data, color + '-',
                         label=label)
        axs[1].plot(tr.times("matplotlib")[imin], tr.data[imin], color='k',
                    marker='x')
        axs[1].plot(tr.times("matplotlib")[imax], tr.data[imax], color='k',
                    marker='x')
        axs[2].plot(tr.times("matplotlib")[imin], tr.data[imin], color='k',
                    marker='x')
        axs[2].plot(tr.times("matplotlib")[imax], tr.data[imax], color='k',
                    marker='x')
        axs[2].set_ylabel(trans_units)

        plt.tight_layout()
        plt.draw()
        plt.show()
        # plt.show(block=False)
        plt.pause(0.001)

    @staticmethod
    def _choose_color(trace, Zcomps='Z3', Ncomps='N1Y', Ecomps='E2X'):
        """
        :param Zcomps: component characters to interpret as "Z" channel
        :param Ncomps: component characters to interpret as "N" channel
        :param Ecomps: component characters to interpret as "E" channel
        """
        comp = trace.stats.channel[-1]
        if comp in Ncomps:
            return 'r', 'N'
        elif comp in Ecomps:
            return 'b', 'E'
        elif comp in Zcomps:
            return 'k', 'Z'
        else:
            return 'g', comp


def pk2pk(stream, start_time, end_time):
    """
    Return the time, period and amplitude of the maximum peak-peak

    :param stream: waveform data
    :param pick_time: time of reference (usually a pick)
    :param start_time: start of amplitude window (UTCDateTime)
    :param end_time: end of amplitude window (UTCDateTime)
    :returns: Amp object
    """
    amp = Amp(value=0)
    if start_time > end_time:
        sta = list(set([tr.stats.station for tr in stream]))
        log(f'station {sta} pk2pk amplitude window start_time > end_time',
            'error')
        return None
    for tr in stream:
        window = tr.copy().trim(start_time, end_time, pad=True).data
        sr = tr.stats.sampling_rate
        # print(f'tr = {tr}')
        # print(f'window = {window}')
        i_mins, _ = find_peaks(-1. * window)
        i_maxs, _ = find_peaks(window)

        if len(i_mins) == 0 or len(i_maxs) == 0:
            return(None)

        # Make sure there is an i_max after the last i_min
        if i_mins[-1] == len(window)-1:
            i_mins = i_mins[:-1]
        if i_maxs[-1] <= i_mins[-1]:
            i_maxs = np.append(i_maxs, len(window)-1)

        for imin in i_mins:
            imax = np.min(i_maxs[i_maxs > imin])
            # print(imin, imax)
            amplitude = window[imax] - window[imin]
            if amplitude > amp.value:
                amp = Amp(amplitude, (imax - imin) / sr,
                          start_time + (imin / sr), tr.stats.channel)
    return amp


def get_response(filename, format=None, component=None):
    """
    Read response file and output PoleZeros object

    :param filename: name of the file to read
    :param format: 'GSE', 'JSON_PZ, 'SACPZ', 'STATIONXML' or None
    :param component: component to read, if STATIONXML
    """
    if format.upper() == 'GSE':
        paz = PAZ.read_gse_response(filename)
    elif format.upper() in ['JSON', 'JSON_PZ']:
        paz = PAZ.read_json_pz(filename)
    elif format.upper() == 'SACPZ':
        paz = PAZ.read_sac_pz(filename)
    elif format.upper() == 'STATIONXML':
        assert component is not None
        paz = PAZ.read_stationxml(filename, '*' + component)
    else:
        paz = PAZ.read_baillard_pz(filename)
    return paz


# def _read_JSON_PZ(filename):
#     """
#     read JSON paz format
#     """
#     with open(filename) as f:
#         paz = json.read(f)
#     A0 = 1
#     for p in paz.poles:
#         p = np.complex(p[0], p[1])
#         A0 *= (2*np.pi*paz['f_ref'] - p)
#     for z in paz.zeros:
#         z = np.complex(z[0], z[1])
#         A0 /= (2*np.pi*paz['f_ref'] - z)
#     A0 = _calc_A0(paz['poles'], paz['zeros'], paz['f_ref'],
#                   paz.get('A0', None))
#     paz['A0'] = A0
#     return paz
#
#
# def _calc_A0(poles, zeros, f_ref, proposed_A0=None):
#     """ Calculation normalization factor """
#     A0 = 1
#     for p in poles:
#         A0 *= (2*np.pi*f_ref - p)
#     for z in zeros:
#         A0 /= (2*np.pi*f_ref - z)
#     A0 = abs(A0)
#     if proposed_A0 is not None:
#         if A0/proposed_A0 < 0.99 or A0/proposed_A0 > 1.01:
#             log('file said A0 at {:.3g} Hz = {:.3g}, calculated {:.3g}'
#                 .format(f_ref, proposed_A0, A0), 'warning')
#     return A0
#
#
# def _read_Baillard_PZ(filename):
#     """
#     read Baillard PZ format
#     """
#     info = {}
#     with open(filename, 'r') as foc:
#         i = 0
#         a, info['poles'], info['zeros'] = [], [], []
#         for tline in foc:
#             if i < 4:
#                 a.append(float(tline))
#             elif i == 4:
#                 n_poles = float(tline)
#             elif i == 5:
#                 n_zeros = float(tline)
#             elif i >= 6 and i < (6 + n_poles):
#                 A = tline.split()
#                 info['poles'].append(np.complex(float(A[0]), float(A[1])))
#             elif i >= 6 + n_poles and i < 6 + n_poles + n_zeros:
#                 A = tline.split()
#                 info['zeros'].append(np.complex(float(A[0]), float(A[1])))
#             else:
#                 info['amplifier_gain'] = float(tline)
#             i += 1
#     info['sensor_gain'] = a[0]
#     info['f_ref'] = 1 / a[1]
#     info['sampling_rate'] = a[2]
#     info['A0'] = a[3]
#     _calc_A0(info['poles'], info['zeros'], info['f_ref'], info['A0'])
#     return info
#
#
# # read_GSE.m
# def _read_GSE(name=None):
#     """
#     read a GSE response file
#
#     :param name:    Name of GSE2 response file
#     :returns: Structure with all parameters
#     """
# # read_GSE.m:6
#     info = dict()
#     with open(name, 'r') as fic:
#         for line in fic:
#             if 'CAL' == line[:3]:
#                 info['station'] = line[6:10]
#                 info['type'] = line[20:29].trim()
#                 info['sensor_gain'] = float(line[29:40])
#                 info['sampling_rate'] = float(line[46:56])
#             else:
#                 if 'PAZ' == line[:3]:
#                     info['A0'] = float(line[11:25])
#                     num_poles = float(line[41:43])
#                     num_zeros = float(line[45])
#                     info['poles'] = []
#                     info['zeros'] = []
#                     for p in np.arange(num_poles):
#                         A = next(fic).split()
#                         info['poles'].extend(np.complex(float(A[0]),
#                                                         float(A[1])))
#                     for z in np.arange(num_zeros):
#                         A = next(fic).split()
#                         info['zeros'].extend(np.complex(float(A[0]),
#                                                         float(A[1])))
#                 else:
#                     if 'DIG' == line[:3]:
#                         info['amplifier_gain'] = float(line[9:23])
#                         info['digitizer'] = line[36:42].trim()
#                     else:
#                         continue
#     info['f_ref'] = 1
#     _calc_A0(info['poles'], info['zeros'], info['f_ref'], info['A0'])
#     return info
