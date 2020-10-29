#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test ps_picker
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA @UnusedWildImport

import os
# import glob
import unittest
import inspect
import difflib
import pprint

from obspy import read as obspy_read
from obspy.core import UTCDateTime
from obspy.core.event.origin import Pick
from obspy.core.event.base import WaveformStreamID

# from obsinfo.misc.info_files import _read_json_yaml
from ps_picker.local_amplitude import (LocalAmplitude, get_response)
from ps_picker.paz import PAZ

pp = pprint.PrettyPrinter(indent=4)


class TestADDONSMethods(unittest.TestCase):
    """
    Test suite for obsinfo operations.
    """
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(inspect.getfile(
            inspect.currentframe())))
        self.testing_path = os.path.join(self.path, "data")

    def assertTextFilesEqual(self, first, second, msg=None):
        with open(first) as f:
            str_a = f.read()
        with open(second) as f:
            str_b = f.read()

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

    def test_paz(self):
        """
        Test calculating amplitudes
        """
        obj = PAZ(poles=[6.228 + 6.228j, 6.228 - 6.228j], zeros=[0, 0],
                     input_units='m/s', output_units='counts',
                     passband_gain=1500., ref_freq=1.)
        obj_nm = PAZ(poles=[6.228 + 6.228j, 6.228 - 6.228j], zeros=[0, 0, 0],
                     input_units='nm', output_units='counts',
                     passband_gain=9.42e-6, ref_freq=1.)
        self.assertNotEqual(obj, obj_nm)
        self.assertAlmostEqual(obj.norm_factor, 2.205, places=3)
        obj.input_units='nm'
        self.assertAlmostEqual(obj.norm_factor, 0.3509, places=3)
        self.assertEqual(obj, obj_nm)

    def test_norm_factor(self):
        """
        Test calculating norm_factor and fixing norm_factor
        """
        obja = PAZ(poles=[6.228 + 6.228j, 6.228 - 6.228j], zeros=[0, 0],
                     input_units='m/s', output_units='counts',
                     passband_gain=1500., ref_freq=1.)
        objb = obja.copy()
        objb.norm_factor = obja.norm_factor
        self.assertAlmostEqual(obja.norm_factor, 2.205, places=3)
        self.assertAlmostEqual(objb.norm_factor, 2.205, places=3)
        obja.input_units = 'nm'
        objb.input_units = 'nm'
        self.assertAlmostEqual(obja.norm_factor, 0.3509, places=3)
        self.assertAlmostEqual(objb.norm_factor, 2.205, places=3)

    def test_resp_file_read(self):
        """
        Test reading response files
        """
        filea = os.path.join(self.testing_path, "SPOBS2_response.txt")
        paza = get_response(filea, '')
        fileb = os.path.join(self.testing_path, "SPOBS2_response.json")
        pazb = get_response(fileb, 'JSON_PZ')
        filec = os.path.join(self.testing_path, "SPOBS2_response.GSE")
        pazc = get_response(filec, 'GSE')
        filed = os.path.join(self.testing_path, "SPOBS2_response.SACPZ")
        pazd = get_response(filed, 'SACPZ')
        pazd.ref_freq = 10.
        pazd.input_units = 'nm'
        filee = os.path.join(self.testing_path, "1T.MOSE.STATION.xml")
        paze = get_response(filee, 'STATIONXML', component='3')
        paze.input_units = 'nm'
        for paz in [pazb, pazc, pazd, paze]:
            self.assertEqual(paza, paz)

    def test_amplitude(self):
        """
        Test calculating amplitudes
        """
        datafile = os.path.join(self.testing_path,
                                "20190519T060917_MONA.mseed")
        respfile = os.path.join(self.testing_path, "SPOBS2_response.json")
        stream = obspy_read(datafile, 'MSEED')
        wid = WaveformStreamID(network_code=stream[0].stats.network,
                               station_code=stream[0].stats.station,
                               channel_code=stream[0].stats.channel)
        Ppick = Pick(time=UTCDateTime('2019-05-19T06:09:48.83'),
                     phase_hint='P', waveform_id=wid)
        Spick = Pick(time=UTCDateTime('2019-05-19T06:09:51.52'),
                     phase_hint='S', waveform_id=wid)
        obj = LocalAmplitude(stream, [Ppick, Spick], respfile, 'JSON_PZ')
        amp_wood_calc = obj.get_iaml(plot=False, method='wood_calc')
        amp_wood_est = obj.get_iaml(plot=False, method='wood_est')
        amp_raw = obj.get_iaml(plot=False, method='raw_disp')
        self.assertAlmostEqual(amp_wood_calc.generic_amplitude,
                               3.7848367595677965e-07)
        self.assertAlmostEqual(amp_wood_calc.period, 0.12)


def suite():
    return unittest.makeSuite(TestADDONSMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
