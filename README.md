# SPSI

## Description:

SPSI (simplest possible structure interpolator) is a tool for creating intermediate macromolecular structures from a pair of structures - otherwise known as morphs.  These can be useful for illustrating conformational changes in macromolecules, and are often displayed as movies.  The startpoint and endpoint models are input as PDB (Protein Data Bank) files, and the output is a multi-model PDB file containing both input models and the interpolated intermediate models.

There are other ways to do this (see the Alternatives section below), which may have more features.  SPSI is an attempt to handle things as simply as possible (in an effort to reduce potential installation and configuration issues).


## Installation:
SPSI requires [python](http://www.python.org)>= 2.5 (earlier versions with sqlite3 may function, but are untested).  Python 3 also works.  SPSI can be used stand-alone (from the command line), or from within [PyMOL](http://pymol.org/).  SPSI has been tested on linux and OS X, but should run on any system with a standard python installation.

### Stand-Alone Installation: 

Copy spsi.py to a directory in your $PATH.

### PyMOL-Integrated Installation:

Copy spsi.py to a directory in PyMOL's $PYTHONPATH (current working directory, or something in sys.path).

## Usage:
For endpoint models A.pdb and B.pdb with 15 interpolated models generated:

### Command line usage: spsi.py A.pdb B.pdb 15

### From within PyMOL:
1. run spsi.py; 
2. load A.pdb
3. load B.pdb
4. morph('A','B','output_morph', nsteps=15) 

The nsteps argument is optional; the default is 10 steps of interpolation.

## Contributers:
Just me (Pete Meyer) so far.

## Alternatives:
RigiMOL (bundled with incentive PyMOL) and LSQMAN (from [Uppsala Software Factory](http://xray.bmc.uu.se/usf/)) are some other (probably more feature rich) options for creating morphs.  There may be others I'm not aware of.

