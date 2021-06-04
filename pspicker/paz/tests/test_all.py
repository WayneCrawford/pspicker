#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test pspicker
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
from pathlib import Path

from obspy import read as obspy_read
from obspy.core import UTCDateTime
from obspy.core.event.origin import Pick, Arrival
from obspy.core.event.magnitude import Amplitude
from obspy.core.event.base import WaveformStreamID, QuantityError
from obspy.core.inventory import read_inventory

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
        self.path = os.path.dirname(os.path.abspath(inspect.getfile(
            inspect.currentframe())))
        self.testing_path = os.path.join(self.path, "data")

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
            # print(str_a)
            # print(str_b)
            first_lines = str_a.splitlines(True)
            second_lines = str_b.splitlines(True)
            delta = difflib.unified_diff(
                first_lines, second_lines,
                fromfile=str(first), tofile=str(second))
            # result = list(delta)
            # print(result)
            message = ''.join(delta)

            if msg:
                message += " : " + msg

            self.fail("Multi-line strings are unequal:\n" + message)

    def test_readwrite(self):
        """
        Test calculating amplitudes
        """
        fmin = 0.001
        inv = read_inventory(str(Path(self.testing_path) / "stations_mayobs.xml"),
                             'STATIONXML')
        for net in inv:
            for sta in net:
                # print(f'{sta.code=}')
                ch = sta.select(channel='*[Z3]')[0]
                info_str = f'{sta.code}.{ch.code}'
                resp = ch.response
                # fig = resp.plot(min_freq=fmin, label=info_str, show=False)
                paz = PAZ.from_obspy_response(resp)
                # paz.plot(min_freq=fmin, label='PAZ', axes=fig.axes, sym='g.',
                #          show=False)
                fname = f'{info_str}.json'
                paz.write_json_pz(fname)
                paz2 = PAZ.read_json_pz(fname)
                paz2.input_units = paz.input_units
                self.assertEqual(paz, paz2)
                self.assertTextFilesEqual(fname, Path(self.testing_path) / fname)
                os.remove(fname)
                # paz2.plot(min_freq=fmin, label='JSON_PAZ', axes=fig.axes,
                #           sym='rx')

    def test_paz(self):
        """
        Test calculating amplitudes
        """
        obj = PAZ.from_refgain(1500., 
                               poles=[6.228 + 6.228j, 6.228 - 6.228j],
                               zeros=[0, 0],
                               ref_freq=1.,
                               input_units='m/s',
                               output_units='counts')
        obj_nm = PAZ.from_refgain(9.42e-6,
                                  poles=[6.228 + 6.228j, 6.228 - 6.228j],
                                  zeros=[0, 0, 0],
                                  ref_freq=1.,
                                  input_units='nm',
                                  output_units='counts')
        self.assertNotEqual(obj, obj_nm)
        self.assertAlmostEqual(obj.norm_factor, 2.205, places=3)
        obj.input_units='nm'
        self.assertAlmostEqual(obj.calc_norm_factor(), 0.3509, places=3)
        self.assertEqual(obj, obj_nm)

    def test_norm_factor(self):
        """
        Test calculating norm_factor
        """
        obja = PAZ.from_refgain(1500.,
                                poles=[6.228 + 6.228j, 6.228 - 6.228j],
                                zeros=[0, 0],
                                ref_freq=1.,
                                input_units='m/s',
                                output_units='counts')
        self.assertAlmostEqual(obja.norm_factor, 2.205, places=3)
        obja.input_units = 'nm'
        self.assertAlmostEqual(obja.calc_norm_factor(), 0.3509, places=3)

    def test_resp_file_read(self):
        """
        Test reading response files
        """
        file_ref = os.path.join(self.testing_path, "SPOBS2_response.txt")
        paz_ref = get_response(file_ref, '')
        suffixes = [".json", ".GSE", '.SACPZ']
        formats = ['JSON_PZ', 'GSE', 'SACPZ']
        for suff, form in zip(suffixes, formats):
            fname = os.path.join(self.testing_path, "SPOBS2_response" + suff)
            paz = get_response(fname, form)
            paz.ref_freq = paz_ref.ref_freq
            paz.input_units = paz_ref.input_units
            self.assertEqual(paz_ref, paz,
                             msg=f'Baillard PAZ != {form}: {paz_ref}, {paz}')
        fname = os.path.join(self.testing_path, "1T.MOSE.STATION.xml")
        paz = get_response(fname, 'STATIONXML', component='3')
        paz.input_units = paz_ref.input_units
        paz.ref_freq = paz_ref.ref_freq
        self.assertEqual(paz_ref, paz,
                         msg=f'Baillard PAZ != STATIONXML: {paz_ref}, {paz}')


def suite():
    return unittest.makeSuite(TestADDONSMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
