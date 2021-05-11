#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test pspicker
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA @UnusedWildImport

import unittest
import inspect
import difflib
import pprint
from pathlib import Path

from obspy import read as obspy_read
from obspy.core import UTCDateTime
from obspy.core.event.origin import Pick, Arrival
from obspy.core.event.magnitude import Amplitude
from obspy.core.event.base import WaveformStreamID, QuantityError

# from obsinfo.misc.info_files import _read_json_yaml
from pspicker.pspicker import PSPicker
from pspicker.local_amplitude import (LocalAmplitude, get_response)
from pspicker.paz import PAZ
from pspicker.logger import setup_log

pp = pprint.PrettyPrinter(indent=4)
setup_log()


class TestADDONSMethods(unittest.TestCase):
    """
    Test suite for obsinfo operations.
    """
    def setUp(self):
        self.path = Path("").resolve().parent
        self.testing_path = Path(self.path) / "tests" / "data"
#         self.path = os.path.dirname(os.path.abspath(inspect.getfile(
#             inspect.currentframe())))
#         self.testing_path = os.path.join(self.path, "data")

    def assertTextFilesEqual(self, first, second, ignore_lines = [],
                             msg=None):
        with open(first) as f:
            str_a = f.read()
        with open(second) as f:
            str_b = f.read()
        if len(ignore_lines) > 0:
            lines = str_a.splitlines(True)
            lines = [l for l in lines if lines.index(l) not in ignore_lines]
            str_a = ''.join(lines)
            lines = str_b.splitlines(True)
            lines = [l for l in lines if lines.index(l) not in ignore_lines]
            str_b = ''.join(lines)
        
        if str_a != str_b:
            first_lines = str_a.splitlines(True)
            second_lines = str_b.splitlines(True)
            delta = difflib.unified_diff(
                first_lines, second_lines,
                fromfile=first, tofile=second)
            message = ''.join(delta)

            if msg:
                message += " : " + msg

            self.fail("Multi-line strings are unequal:\n" + message)

    def test_resp_file_read(self):
        """
        Test reading response files
        """
        filea = str(self.testing_path / "SPOBS2_response.txt")
        paza = get_response(filea, '')
        filec = str(self.testing_path / "SPOBS2_response.GSE")
        pazc = get_response(filec, 'GSE')
        filed = str(self.testing_path / "SPOBS2_response.SACPZ")
        pazd = get_response(filed, 'SACPZ')
        pazd.ref_freq = 10.
        pazd.input_units = 'nm'
        filee = str(self.testing_path / "1T.MOSE.STATION.xml")
        paze = get_response(filee, 'STATIONXML', component='3')
        paze.input_units = 'nm'
        for paz in [pazc, pazd, paze]:
            self.assertEqual(paza, paz)

    def test_amplitude(self):
        """
        Test calculating amplitudes
        """
        datafile = str(self.testing_path / "20190519T060917_MONA.mseed")
        respfile = str(self.testing_path / "SPOBS2_response.json")
        stream = obspy_read(datafile, 'MSEED')
        wid = WaveformStreamID(network_code=stream[0].stats.network,
                               station_code=stream[0].stats.station,
                               channel_code=stream[0].stats.channel)
        Ppick = Pick(time=UTCDateTime('2019-05-19T06:09:48.83'),
                     phase_hint='P', waveform_id=wid)
        Spick = Pick(time=UTCDateTime('2019-05-19T06:09:51.52'),
                     phase_hint='S', waveform_id=wid)
        obj = LocalAmplitude(stream, [Ppick, Spick], respfile, 'JSON_PZ')
        amp_wood_calc, _ = obj.get_iaml(plot=False, method='wood_calc')
        amp_wood_est, _ = obj.get_iaml(plot=False, method='wood_est')
        amp_raw, _ = obj.get_iaml(plot=False, method='raw_disp')
        self.assertAlmostEqual(amp_wood_calc.generic_amplitude,
                               1097.555091056036)
        self.assertAlmostEqual(amp_wood_calc.period, 0.112)


    def test_nordic_write(self):
        """
        Test calculating amplitudes
        """
        compare_file = str(self.testing_path / "test.nordic")
        otime = UTCDateTime('2019-05-19T06:09:48')
        wid = WaveformStreamID(network_code='4G', station_code='STAT')
        wid.channel_code = 'SHZ'
        Ppick = Pick(time=otime + 2.2, phase_hint='P', waveform_id=wid,
                      evaluation_mode='automatic',
                      time_errors=QuantityError(0.01))
        Parrival = Arrival(pick_id=Ppick.resource_id, time_weight=0)
        wid.channel_code = 'SH1'
        Spick = Pick(time=otime+3.5, phase_hint='S', waveform_id=wid,
                      evaluation_mode='automatic',
                      time_errors=QuantityError(0.05))
        wid.channel_code = 'SH2'
        Sarrival = Arrival(pick_id=Spick.resource_id, time_weight=1)
        Apick = Pick(time=otime + 4.0, phase_hint='IAML', waveform_id=wid,
                      evaluation_mode='automatic',
                      time_errors=QuantityError(0.1))
        Amp = Amplitude(generic_amplitude=410.e-9, type='IAML', unit='m',
                        period=0.22, magnitude_hint='ML', category='period',
                        pick_id=Apick.resource_id,
                        waveform_id=Apick.waveform_id)
        picks = [Ppick, Spick, Apick]
        amps = [Amp]
        arrivals = [Parrival, Sarrival]
        PSPicker.save_nordic_event(picks, otime, '.', 'test.nordic',
                                   amps, arrivals=arrivals, debug=False)
        self.assertTextFilesEqual('test.nordic', compare_file, ignore_lines=[1])
        Path("test.nordic").unlink()
        for p in Path(".").glob('run_*.log'):
            p.unlink()

def suite():
    return unittest.makeSuite(TestADDONSMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
