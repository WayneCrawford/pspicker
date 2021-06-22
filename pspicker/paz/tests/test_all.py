#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test PAZ
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

from obspy.core.inventory import read_inventory

from pspicker.local_amplitude import get_response
from pspicker.paz import PAZ
# from pspicker.logger import setup_log

pp = pprint.PrettyPrinter(indent=4)
# setup_log()


class TestADDONSMethods(unittest.TestCase):
    """
    Test suite for obsinfo operations.
    """
    def setUp(self):
        self.path = Path(inspect.getfile(inspect.currentframe())).resolve().parent
        self.data_path = Path(self.path) / "data"

    def assertTextFilesEqual(self, first, second, ignore_lines=[], msg=None):
        with open(first) as f:
            str_a = f.read()
        with open(second) as f:
            str_b = f.read()
        if len(ignore_lines) > 0:
            lines = str_a.splitlines(True)
            lines = [x for x in lines if lines.index(x) not in ignore_lines]
            str_a = ''.join(lines)
            lines = str_b.splitlines(True)
            lines = [x for x in lines if lines.index(x) not in ignore_lines]
            str_b = ''.join(lines)
        if str_a != str_b:
            first_lines = str_a.splitlines(True)
            second_lines = str_b.splitlines(True)
            delta = difflib.unified_diff(
                first_lines, second_lines,
                fromfile=str(first), tofile=str(second))
            message = ''.join(delta)
            if msg:
                message += " : " + msg
            self.fail("Multi-line strings are unequal:\n" + message)

    def test_readwrite(self):
        """
        Test calculating amplitudes
        """
        # fmin = 0.001
        inv = read_inventory(str(self.data_path / "stations_mayobs.xml"),
                             'STATIONXML')
        for net in inv:
            for sta in net:
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
                self.assertTextFilesEqual(fname, self.data_path / fname)
                os.remove(fname)
                # paz2.plot(min_freq=fmin, label='JSON_PAZ', axes=fig.axes,
                #           sym='rx')

    def test_ampls(self):
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
        obj.input_units = 'nm'
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

    def test_file_reads(self):
        """
        Test reading response files
        """
        paz=[]
        paz.append(PAZ.read_baillard_pz(self.data_path / "SPOBS2_response.txt"))
        paz.append(PAZ.read_json_pz(self.data_path / "SPOBS2_response.json"))
        paz.append(PAZ.read_gse_response(self.data_path / "SPOBS2_response.GSE"))
        paz.append(PAZ.read_sac_pz(self.data_path / "SPOBS2_response.SACPZ"))
        paz[-1].ref_freq = 10.
        paz.append(PAZ.read_stationxml(str(self.data_path / "1T.MOSE.STATION.xml"),
                                       '*3'))
        for pz in paz[1:]:
            pz.input_units = 'nm'
            self.assertEqual(paz[0], pz)

    def test_sacpz_readwrite(self):
        """
        Test reading and writing to sac pz file
        """
        in_file = self.data_path / "SPOBS2_response.SACPZ"
        out_file = Path("SPOBS2_response.SACPZ")
        paz = PAZ.read_sac_pz(in_file)
        paz.write_sac_pz(out_file)
        self.assertTextFilesEqual(in_file, out_file)
        out_file.unlink()


def suite():
    return unittest.makeSuite(TestADDONSMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
