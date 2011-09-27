#!/usr/bin/env python

import unittest

from spsi import parse_pdb_atom_line, linear_interpolate, generate_pdb_atom_line

class TestIO(unittest.TestCase):
    def test1(self):
        ''' atom parse 1 '''
        ln = 'ATOM      1  N   VAL A   2     -10.675 -27.377 -28.169  1.00 93.06      A    N\n'
        d = parse_pdb_atom_line(ln)
        self.assert_( d['serial'] == 1)
        self.assert_( d['x'] == -10.675)
        self.assert_( d['y'] == -27.377)
        self.assert_( d['z'] == -28.169)
        self.assert_( d['chainID'] == 'A')
        self.assert_( d['name'] == ' N  ')
        self.assert_( d['resSeq'] == 2)
    def test2(self):
        ''' atom parse 2 '''
        ln = 'ATOM  18027  CA  VAL I  43      32.664-109.147  -0.440  1.00139.76      I    C'
        d = parse_pdb_atom_line(ln)
        self.assert_( d['x'] == 32.664)
        self.assert_( d['y'] == -109.147)
        self.assert_( d['z'] == -0.440)
        self.assert_( d['name'] == ' CA ')
        self.assert_( d['resSeq'] == 43)
        self.assert_( d['chainID'] == 'I')
    def test2b(self):
        ''' atom parse 3 '''
        ln = 'ATOM   7838  N   LEU A1000     -31.768 -60.286 -16.708  1.00103.83      A    N'
        d = parse_pdb_atom_line(ln)
        self.assert_( d['x'] == -31.768)
        self.assert_( d['y'] == -60.286)
        self.assert_( d['z'] == -16.708)
        self.assert_( d['name'] == ' N  ')
        self.assert_( d['resSeq'] == 1000)
        self.assert_( d['chainID'] == 'A')
    def test3(self):
        ''' interpolation test '''
        pcts = [ 0.0, 0.25, 0.5, 0.75, 1.0 ]
        st = -5.0
        en = 5.0
        pts = [-5.0, -2.5, 0.0, 2.5, 5.0 ]
        ipts = [ linear_interpolate(st, en, p) for p in pcts ]
        self.assert_( pts == ipts )
    def test4(self):
        ''' output formatting test 1'''
        ln1 = 'ATOM  18027  CA  VAL I  43      32.664-109.147  -0.440  1.00139.76      I    C\n'
        ln2 = 'ATOM      1  CA  VAL I  43      32.664-109.147  -0.440  1.00 23.42\n'
        a = parse_pdb_atom_line(ln1)
        oln = generate_pdb_atom_line(1, a['name'], a['resName'], a['chainID'], a['resSeq'], a['x'], a['y'], a['z'] )
        self.assert_( oln == ln2 )
    def test5(self):
        ''' output formatting test 2'''
        ln1 = 'ATOM   7838  N   LEU A1000     -31.768 -60.286 -16.708  1.00103.83      A    N\n'
        ln2 = 'ATOM      1  N   LEU A1000     -31.768 -60.286 -16.708  1.00 23.42\n'
        a = parse_pdb_atom_line(ln1)
        oln = generate_pdb_atom_line(1, a['name'], a['resName'], a['chainID'], a['resSeq'], a['x'], a['y'], a['z'] )
        self.assert_( oln == ln2 )
    def test6(self):
        ''' output formatting test 3'''
        ln1 = 'ATOM   7016  CB  SER A1449     -35.964 -34.300 -27.433  1.00100.34      A    C'
        ln2 = 'ATOM      1  CB  SER A1449     -35.964 -34.300 -27.433  1.00 23.42\n'
        a = parse_pdb_atom_line(ln1)
        oln = generate_pdb_atom_line(1, a['name'], a['resName'], a['chainID'], a['resSeq'], a['x'], a['y'], a['z'] )
        self.assert_( oln == ln2 )

if __name__ == '__main__':
    def getsuite(tc):
        ''' wrapper method to load a test suite '''
        return unittest.TestLoader().loadTestsFromTestCase(tc)
    suite = getsuite(TestIO)
    unittest.TextTestRunner(verbosity=2).run(suite)
